classdef DEIF_FormationEstimationByRangeAngleWithReference < ...
    DIF_LinearDynamics & ...
    MultiagentUtilityBase

    properties (SetAccess = private)
        range_sensor_ % Instance of RangeMeasurementMultiAgentWithReference class
        angle_sensor_ % Instance of AngleMeasurementMultiAgentWithReference class
        counter
    end
    properties (SetAccess = protected)
        state_vector
        state_covmat
    end
    properties (SetAccess = immutable)
        % Abstract properties of ExtendedInformationFilter
        process_noise_covmat
        discrete_system_matrix
        % Abstract properties of MultiagentUtilityBase
        num_agents
        num_dimensions
        % 
        wait_steps
    end

    methods (Access = public)
        function obj = DEIF_FormationEstimationByRangeAngleWithReference(args)
            obj@DIF_LinearDynamics(args);
            obj.range_sensor_           = RangeMeasurementMultiAgentWithReference(args.range_sensor);
            obj.angle_sensor_           = AngleMeasurementMultiAgentWithReference(args.angle_sensor);
            obj.num_agents              = args.num_agents;
            obj.num_dimensions          = args.num_dimensions;
            obj.state_vector            = args.state_vector;
            obj.state_covmat            = obj.createStateCovarianceMatrix(args.sigma_position, args.sigma_velocity);
            obj.process_noise_covmat    = args.process_noise_covmat;
            obj.discrete_system_matrix  = obj.createDiscreteSystemMatrix(args.discrete_system_matrix);
            obj.counter                 = 0;
            obj.wait_steps              = args.wait_steps;
        end
    end

    methods (Access = public)
        function executeFiltering(this, measures, adjacent_matrix, ...
            outsource_info_vector, outsource_info_matrix, position_ref)
            this.predictStateVectorAndCovariance();
            this.resetObservationInformation();
            this.processMeasurements(measures, adjacent_matrix, position_ref);
            this.clearOutSourceInformationPool();
            this.addOutSourceInformationIntoPool(...
                outsource_info_vector, outsource_info_matrix);
            this.integrateMultiSourceInformation();
            if (this.counter >= this.wait_steps)
                this.computePosteriorPdf();
                % this.clearOutSourceInformationPool();
                this.counter = 0;
            else
                this.counter = this.counter + 1;
            end
        end
    end

    methods (Access = protected)
        function processMeasurements(this, measures, adjacent_matrix, position_ref)
            % measures should have 'ranges' and 'angles' field
            positions = this.getPositionVector();

            % Range measurements
            this.range_sensor_.computeMeasurementVector(positions, position_ref, false);
            this.range_sensor_.setObservationMatrix(positions, position_ref);
            this.range_sensor_.updateMeasurementCovarianceMatrix(adjacent_matrix.range);
            obs_matrix_range = this.range_sensor_.getObservationMatrix();
            obs_covmat_range = this.range_sensor_.getMeasureCovarinaceMatrix();
            measures_predicted_range = this.range_sensor_.getMeasurements();
            this.addObservationInformation(...
                obs_matrix_range, obs_covmat_range, measures.ranges, measures_predicted_range);
            
            % Angle measurements
            this.angle_sensor_.computeMeasurementVector(positions, position_ref, false);
            this.angle_sensor_.setObservationMatrix(positions, position_ref);
            this.angle_sensor_.setMeasurementCovarianceMatrix(adjacent_matrix.angle);
            obs_matrix_angle = this.angle_sensor_.getObservationMatrix();
            obs_covmat_angle = this.angle_sensor_.getMeasureCovarinaceMatrix();
            measures_predicted_angle = this.angle_sensor_.getMeasurements();
            % TODO: How to tuckle the following singular point problem?
            diff_measures = measures_predicted_angle - measures.angles;
            for iMeasures = 1:length(diff_measures)
                if (abs(diff_measures(iMeasures,1)) >= pi)
                    if (measures_predicted_angle(iMeasures,1) > measures.angles(iMeasures,1))
                        measures_predicted_angle(iMeasures,1) = measures_predicted_angle(iMeasures,1) - 2*pi;
                    else
                        measures_predicted_angle(iMeasures,1) = measures_predicted_angle(iMeasures,1) + 2*pi;
                    end
                end
            end
            this.addObservationInformation(...
                obs_matrix_angle, obs_covmat_angle, measures.angles, measures_predicted_angle);   
        end
    end

    methods (Access = protected)
        function addObservationInformation(this,...
            obs_matrix, obs_covmat, measures, measures_predicted)
            H = obs_matrix;
            R = obs_covmat;
            y = measures;
            y_hat = measures_predicted;
            x_hat = this.state_vector;
            this.obs_info_vector = this.obs_info_vector + H.'/R*(y - y_hat + H*x_hat);
            this.obs_info_matrix = this.obs_info_matrix + H.'/R*H;
        end
    end

    methods (Access = private)
        % Getters -------------------------------------------
        function output = getPositionVector(this)
            num_dims = this.num_dimensions;
            output = zeros(num_dims*this.num_agents,1);
            for iAgents = 1:this.num_agents
                output(num_dims*(iAgents-1)+1:num_dims*iAgents,1) ...
                    = this.state_vector(2*num_dims*(iAgents-1)+1:2*num_dims*(iAgents-1)+num_dims,1);
            end
        end
    end

end