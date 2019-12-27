classdef EIF_FormationEstimationByRange < ...
    ExtendedInformationFilter & ...
    RangeMeasurementInterAgents
    properties (SetAccess = private)
    end

    methods
        function obj = EIF_FormationEstimationByRange(args)
            obj@ExtendedInformationFilter(args.eif);
            obj@RangeMeasurementInterAgents(args.rmia);
        end

        function executeInformationFilter(this, measures, adjacent_matrix)
            this.predictStateVectorAndCovariance();
            this.convertMomentsToInformationForm();

            % Range measurements
            positions = this.getPositionVector();
            this.setMeasurementVectorWithoutNoise(positions);
            this.setObservationMatrix(positions);
            this.updateMeasurementCovarianceMatrix(adjacent_matrix);
            obs_matrix = this.getObservationMatrix();
            obs_covmat = this.getMeasureCovarinaceMatrix();
            measures_predicted = this.getMeasurements();
            this.addObservationInformation(...
                obs_matrix, obs_covmat, measures, measures_predicted);
            
            % Position measurements

            this.convertInformationToMomentsForm();
        end

        % Setters -------------------------------------------
        function setStateCovarianceMatrix(this, args)
            P = zeros(size(this.state_covmat));
            num_dims = this.num_dims;
            for iAgents = 1:this.num_agents
                for iDims = 1:num_dims
                    P(2*num_dims*(iAgents-1)+iDims, 2*num_dims*(iAgents-1)+iDims) = args.position_sigma^2;
                    P(2*num_dims*(iAgents-1)+num_dims+iDims, 2*num_dims*(iAgents-1)+num_dims+iDims) = args.velocity_sigma^2;
                end
            end
            this.state_covmat = P;
        end

        function setDiscreteSystemMatrix(this, discrete_system_matrix)
            num_agents = this.num_agents;
            num_vars = this.num_variables;
            num_dims = this.num_dims;
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