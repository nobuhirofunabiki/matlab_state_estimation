classdef UKF_FormationEstimationByRangeAngleWithReference < ...
    UKF_LinearDynamics & ...
    MultiagentUtilityBase

    properties (SetAccess = protected)
        state_vector
        state_covmat
        discrete_system_matrix
        process_noise_covmat
    end
    properties (SetAccess = private)
        range_sensor_       % Instance of RangeMeasurementMultiAgentWithReference class
        angle_sensor_       % Instance of AngleMeasurementMultiAgentWithReference class
        % Ratio actual noise vs noise model in the estimator
        ratio_range_noise
        ratio_angle_noise
    end
    properties (SetAccess = immutable)
        num_dimensions
        num_agents
    end

    methods (Access = public)
        function obj = UKF_FormationEstimationByRangeAngleWithReference(args)
            obj@UKF_LinearDynamics(args);
            obj.range_sensor_           = RangeMeasurementMultiAgentWithReference(args.range_sensor);
            obj.angle_sensor_           = AngleMeasurementMultiAgentWithReference(args.angle_sensor);
            obj.ratio_range_noise       = args.ratio_range_noise;
            obj.ratio_angle_noise       = args.ratio_angle_noise;
            obj.num_agents              = args.num_agents;
            obj.num_dimensions          = args.num_dimensions;
            obj.state_vector            = args.init_state_vector;
            obj.process_noise_covmat    = args.process_noise_covmat;
            obj.discrete_system_matrix  = obj.createDiscreteSystemMatrix(args.discrete_system_matrix);
            obj.state_covmat            = obj.createStateCovarianceMatrix(args.sigma_position, args.sigma_velocity);
        end
    end

    methods (Access = public)
        function executeUnscentedKalmanFilter(this, measurements, adjacent_matrix, position_ref)
            this.sigma_points = this.generateSigmaPoints();
            this.predictStates();
            this.sigma_points = this.generateSigmaPoints();
            this.processMeasurements(measurements, adjacent_matrix, position_ref);
            this.updateByMeasurements();
        end
    end

    methods (Access = protected)
        function processMeasurements(this, measurements, adjacent_matrix, position_ref)
            this.measurements   = vertcat(measurements.ranges, measurements.angles);
            num_sigma_points    = size(this.sigma_points, 2);
            num_measures        = size(this.measurements, 1);
            num_ranges          = size(measurements.ranges, 1);
            num_angles          = size(measurements.angles, 1);
            z_pred              = zeros(num_measures,1);
            Z                   = zeros(num_measures,num_sigma_points);
            obs_covmat          = zeros(num_measures, num_measures);
            this.range_sensor_.updateMeasurementCovarianceMatrix(adjacent_matrix.range);
            this.angle_sensor_.setMeasurementCovarianceMatrix(adjacent_matrix.angle);
            obs_covmat(1:num_ranges, 1:num_ranges) ...
                = this.range_sensor_.getMeasureCovarinaceMatrix() * this.ratio_range_noise;
            obs_covmat(num_ranges+1:num_ranges+num_angles, num_ranges+1:num_ranges+num_angles) ...
                = this.angle_sensor_.getMeasureCovarinaceMatrix() * this.ratio_angle_noise;
            for iPoints = 1:num_sigma_points
                positions = this.getPositionVector(iPoints);
                this.range_sensor_.computeMeasurementVector(positions, position_ref, false);
                this.angle_sensor_.computeMeasurementVector(positions, position_ref, false);
                measures_predicted_angle = this.angle_sensor_.getMeasurements();
                % TODO: How to tuckle the following singular point problem?
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
                Z(1:num_ranges, iPoints)                            = this.range_sensor_.getMeasurements();
                Z(num_ranges+1:num_ranges+num_angles, iPoints)      = measures_predicted_angle;
                z_pred = z_pred + this.weights_mean(1,iPoints)*Z(:,iPoints);
            end
            this.z_pred = z_pred;
            this.Z_diff = Z - z_pred(:,ones(1,num_sigma_points));
            this.Pzz    = this.Z_diff*diag(this.weights_cov)*(this.Z_diff)' + obs_covmat;
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