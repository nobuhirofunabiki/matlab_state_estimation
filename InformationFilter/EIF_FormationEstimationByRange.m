classdef EIF_FormationEstimationByRange < ...
    ExtendedInformationFilter & ...
    MultiagentUtilityBase

    properties (SetAccess = immutable)
        % Abstract properties of ExtendedInformationFilter
        process_noise_covmat
        discrete_system_matrix
        % Abstract properties of MultiagentUtilityBase
        num_agents
        num_dimensions
    end
    properties (SetAccess = protected)
        % Abstract properties of InformationFilterBase
        state_vector
        state_covmat
    end
    properties (SetAccess = private)
        range_sensor_       % Instance of RangeMeasurementInterAgents class
        angle_sensor_       % Instance of AngleMeasurementMultiAgents class
        position_sensor_    % Instance of PositionMeasurementMultiAgents class
    end

    methods (Access = public)
        function obj = EIF_FormationEstimationByRange(args)
            obj@ExtendedInformationFilter(args);
            obj.checkConstructorArguments(args);
            obj.range_sensor_           = RangeMeasurementInterAgents(args.rmia);
            obj.angle_sensor_           = AngleMeasurementMultiAgents(args.angle);
            obj.position_sensor_        = PositionMeasurementMultiAgents(args.pmb);
            obj.num_agents              = args.num_agents;
            obj.num_dimensions          = args.num_dimensions;
            obj.state_vector            = args.state_vector;
            obj.process_noise_covmat    = args.process_noise_covmat;
            obj.discrete_system_matrix  = obj.createDiscreteSystemMatrix(args.discrete_system_matrix);
            obj.state_covmat            = obj.createStateCovarianceMatrix(args.sigma_position, args.sigma_velocity);
        end
    end

    methods (Access = private)
        function checkConstructorArguments(this, args)
            num_vars = this.num_variables;
            disp("Check constructor arguments for EIF_FormationEstimationByRange");
            assert(isequal(size(args.state_vector), [num_vars, 1]), ...
                "state_vector is NOT a correct size vector");
            assert(isequal(size(args.process_noise_covmat), [num_vars, num_vars]), ...
                "process_noise_covmat is NOT a correct size vector");
        end
    end

    methods (Access = public)
        function executeInformationFilter(this, measures, adjacent_matrix)
            % measures should have 'ranges' and 'positions' field
            this.predictStateVectorAndCovariance();
            this.convertMomentsToInformationForm();

            positions = this.getPositionVector();

            % Range measurements
            this.range_sensor_.calculateMeasurementVectorWithoutNoise(positions);
            this.range_sensor_.setObservationMatrix(positions);
            this.range_sensor_.updateMeasurementCovarianceMatrix(adjacent_matrix.range);
            obs_matrix_range = this.range_sensor_.getObservationMatrix();
            obs_covmat_range = this.range_sensor_.getMeasureCovarinaceMatrix();
            measures_predicted_range = this.range_sensor_.getMeasurements();
            this.addObservationInformation(...
                obs_matrix_range, obs_covmat_range, measures.ranges, measures_predicted_range);
            
            % Angle measurements
            this.angle_sensor_.computeMeasurementVector(positions, false);
            this.angle_sensor_.setObservationMatrix(positions);
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

            % Position measurements
            % TODO: Should I use the linear version for addObservationInformation?
            this.position_sensor_.computeMeasurementVector(positions, false);
            obs_matrix_pos = this.position_sensor_.getObservationMatrix();
            obs_covmat_pos = this.position_sensor_.getMeasureCovarinaceMatrix();
            measures_predicted_pos = this.position_sensor_.getMeasurements();
            this.addObservationInformation(...
                obs_matrix_pos, obs_covmat_pos, measures.positions, measures_predicted_pos);

            this.convertInformationToMomentsForm();
        end

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