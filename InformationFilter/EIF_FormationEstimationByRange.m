classdef EIF_FormationEstimationByRange < ExtendedInformationFilter
    
    properties (SetAccess = private)
        range_              % Instance of RangeMeasurementInterAgents class
        position_sensor_    % Instance of PositionMeasurementMultiAgents class
        num_agents
    end

    methods (Access = public)
        function obj = EIF_FormationEstimationByRange(args)
            obj@ExtendedInformationFilter(args.eif);
            obj.checkConstructorArguments(args);
            obj.range_                  = RangeMeasurementInterAgents(args.rmia);
            obj.position_sensor_        = PositionMeasurementMultiAgents(args.pmb);
            obj.num_agents              = args.num_agents;
            obj.state_vector            = args.state_vector;
            obj.process_noise_covmat    = args.process_noise_covmat;
            obj.setStateCovarianceMatrix(args.sigma_position, args.sigma_velocity);
            obj.setDiscreteSystemMatrix(args.discrete_system_matrix);
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
            this.range_.calculateMeasurementVectorWithoutNoise(positions);
            this.range_.setObservationMatrix(positions);
            this.range_.updateMeasurementCovarianceMatrix(adjacent_matrix);
            obs_matrix_range = this.range_.getObservationMatrix();
            obs_covmat_range = this.range_.getMeasureCovarinaceMatrix();
            measures_predicted_range = this.range_.getMeasurements();
            this.addObservationInformation(...
                obs_matrix_range, obs_covmat_range, measures.ranges, measures_predicted_range);
            
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

        % Setters -------------------------------------------
        function setStateCovarianceMatrix(this, sigma_position, sigma_velocity)
            num_dims = this.num_dimensions;
            num_vars = this.num_variables;
            P = zeros(num_vars, num_vars);
            for iAgents = 1:this.num_agents
                for iDims = 1:num_dims
                    P(2*num_dims*(iAgents-1)+iDims, 2*num_dims*(iAgents-1)+iDims) = sigma_position^2;
                    P(2*num_dims*(iAgents-1)+num_dims+iDims, 2*num_dims*(iAgents-1)+num_dims+iDims) = sigma_velocity^2;
                end
            end
            this.state_covmat = P;
        end

        function setDiscreteSystemMatrix(this, discrete_system_matrix)
            num_agents = this.num_agents;
            num_vars = this.num_variables;
            num_dims = this.num_dimensions;
            Ad = zeros(num_vars, num_vars);
            for iAgents = 1:num_agents
                Ad(1+2*num_dims*(iAgents-1):2*num_dims*iAgents,...
                    1+2*num_dims*(iAgents-1):2*num_dims*iAgents)...
                = discrete_system_matrix;
            end
            this.discrete_system_matrix = Ad;
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