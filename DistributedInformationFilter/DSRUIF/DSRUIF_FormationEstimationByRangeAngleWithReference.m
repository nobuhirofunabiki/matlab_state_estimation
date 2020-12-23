classdef DSRUIF_FormationEstimationByRangeAngleWithReference < ...
    DecentralizedSquareRootUnscentedInformationFilterBase & ...
    MultiagentUtilityBase

    properties (SetAccess = private)
        range_sensor_ % Instance of RangeMeasurementMultiAgentWithReference class
        angle_sensor_ % Instance of AngleMeasurementMultiAgentWithReference class
        counter
    end
    properties (SetAccess = protected)
        state_vector
        state_covmat
        sqrtQ
        sqrt_info_matrix
    end
    properties (SetAccess = immutable)
        % Abstract properties of DecentralizedUnscentedInformationFilterBase
        process_noise_covmat
        discrete_system_matrix
        % Abstract properties of MultiagentUtilityBase
        num_agents
        num_dimensions
        % 
        wait_steps
        % Ratio actual noise vs noise model in the estimator
        ratio_range_noise
        ratio_angle_noise
    end

    methods (Access = public)
        function obj = DSRUIF_FormationEstimationByRangeAngleWithReference(args)
            obj@DecentralizedSquareRootUnscentedInformationFilterBase(args);
            obj.range_sensor_           = RangeMeasurementMultiAgentWithReference(args.range_sensor);
            obj.angle_sensor_           = AngleMeasurementMultiAgentWithReference(args.angle_sensor);
            obj.ratio_range_noise       = args.ratio_range_noise;
            obj.ratio_angle_noise       = args.ratio_angle_noise;
            obj.num_agents              = args.num_agents;
            obj.num_dimensions          = args.num_dimensions;
            obj.state_vector            = args.state_vector;
            obj.state_covmat            = obj.createStateCovarianceMatrix(args.sigma_position, args.sigma_velocity);
            obj.process_noise_covmat    = args.process_noise_covmat;
            obj.discrete_system_matrix  = obj.createDiscreteSystemMatrix(args.discrete_system_matrix);
            obj.counter                 = 0;
            obj.wait_steps              = args.wait_steps;
            sqrt_cov_matrix             = chol(obj.state_covmat); 
            obj.sqrt_info_matrix        = transpose(inv(sqrt_cov_matrix));
            obj.sqrtQ                   = chol(args.process_noise_covmat);
            % obj.info_vector             = obj.sqrt_info_matrix*obj.sqrt_info_matrix.'*obj.state_vector;
            % obj.info_vector             = obj.sqrt_info_matrix.'*obj.sqrt_info_matrix*obj.state_vector;
            obj.info_vector             = inv(obj.state_covmat) * obj.state_vector;
        end
    end

    methods (Access = public)
        function executeFiltering(this, measures, adjacent_matrix, ...
            outsource_info_vector, sqrt_outsource_info_matrix, position_ref)
            this.predictStateVectorAndCovariance();
            this.resetObservationInformation();
            this.processMeasurements(measures, adjacent_matrix, position_ref);
            this.clearOutSourceInformationPool();
            this.addOutSourceInformationIntoPool(...
                outsource_info_vector, sqrt_outsource_info_matrix);
            this.integrateMultiSourceInformation();
            if (this.counter >= this.wait_steps)
                this.computePosteriorPdf();
                this.counter = 0;
            else
                this.counter = this.counter + 1;
            end
        end
    end

    methods (Access = protected)
        function processMeasurements(this, measurements, adjacent_matrix, position_ref)
            sqrt_cov_matrix = transpose(inv(this.sqrt_info_matrix));
            this.sigma_points = this.generateSigmaPoints(this.state_vector, sqrt_cov_matrix);

            agg_measurements    = vertcat(measurements.ranges, measurements.angles);
            num_measures        = size(agg_measurements, 1);
            num_ranges          = size(measurements.ranges, 1);
            num_angles          = size(measurements.angles, 1);
            z_pred              = zeros(num_measures,1);
            obs_covmat          = zeros(num_measures, num_measures);
            this.range_sensor_.updateMeasurementCovarianceMatrix(adjacent_matrix.range);
            this.angle_sensor_.setMeasurementCovarianceMatrix(adjacent_matrix.angle);
            obs_covmat(1:num_ranges, 1:num_ranges) ...
                = this.range_sensor_.getMeasureCovarinaceMatrix() * this.ratio_range_noise;
            obs_covmat(num_ranges+1:num_ranges+num_angles, num_ranges+1:num_ranges+num_angles) ...
                = this.angle_sensor_.getMeasureCovarinaceMatrix() * this.ratio_angle_noise;
            this.sqrtR = chol(obs_covmat);
            
            num_sigma_points = length(this.sigma_points);
            Z = zeros(num_measures,num_sigma_points);

            for iPoints = 1:num_sigma_points
                positions = this.getPositionVector(iPoints);
                this.range_sensor_.computeMeasurementVector(positions, position_ref, false);
                this.angle_sensor_.computeMeasurementVector(positions, position_ref, false);
                measures_predicted_angle = this.angle_sensor_.getMeasurements();
                % TODO: How to tackle the following singular point problem?
                diff_measures = measures_predicted_angle - measurements.angles;
                for iMeasures = 1:length(diff_measures)
                    if (abs(diff_measures(iMeasures,1)) >= pi)
                        if (measures_predicted_angle(iMeasures,1) > measurements.angles(iMeasures,1))
                            measures_predicted_angle(iMeasures,1) = measures_predicted_angle(iMeasures,1) - 2*pi;
                        else
                            measures_predicted_angle(iMeasures,1) = measures_predicted_angle(iMeasures,1) + 2*pi;
                        end
                    end
                end
                Z(1:num_ranges, iPoints) = this.range_sensor_.getMeasurements();
                Z(num_ranges+1:num_ranges+num_angles, iPoints) = measures_predicted_angle;
                z_pred = z_pred + this.weights_mean(1,iPoints)*Z(:,iPoints);
            end
            Z_diff = Z - z_pred(:,ones(1,num_sigma_points));
            this.X_diff = this.sigma_points - this.state_vector;
            % Pxy = this.X_diff*diag(this.weights_cov)*Z_diff.';
            Pxy = this.X_diff*diag(this.weights_cov)*Z_diff';
            this.addObservationInformation(Pxy, agg_measurements, z_pred);
        end

        function addObservationInformation(this, Pxy, measures, measures_predicted)
            % sqrtR = transpose(inv(this.sqrtR));
            sqrtR = this.sqrtR;
            y = measures;
            y_hat = measures_predicted;
            x_hat = this.state_vector;
            Sz = this.sqrt_info_matrix;
            Si = sqrtR.' * Pxy.' * Sz * Sz.';
            obs_info_matrix_Si = Si.' * Si;
            this.obs_info_vector = Si.' * sqrtR * (y - y_hat) + Si.' * Si * x_hat;
            Si_square = this.functions_.qr_decomposition_pure(Si.');
            this.sqrt_obs_info_matrix = Si_square;

            info_matrix = Sz.'*Sz;
            obs_info_matrix = Si_square.' * Si_square;
        end
    end

    methods (Access = private)
        function output = getPositionVector(this, iSigmaPoints)
            num_dims = this.num_dimensions;
            output = zeros(num_dims*this.num_agents,1);
            for iAgents = 1:this.num_agents
                output(num_dims*(iAgents-1)+1:num_dims*iAgents,1) ...
                    = this.sigma_points(2*num_dims*(iAgents-1)+1:2*num_dims*(iAgents-1)+num_dims,iSigmaPoints);
            end
        end
    end
    
end