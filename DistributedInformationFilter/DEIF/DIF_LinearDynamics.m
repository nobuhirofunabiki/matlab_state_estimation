classdef DIF_LinearDynamics < DecentralizedInformationFilterBase
    properties (Abstract = true, SetAccess = immutable)
        num_dimensions
        discrete_system_matrix
        process_noise_covmat
    end

    methods (Access = protected)
        function obj = DIF_LinearDynamics(args)
            obj@DecentralizedInformationFilterBase(args);
            obj.process_noise_covmat    = args.process_noise_covmat;
        end
    end

    methods (Access = private)
        function checkConstructorArguments(this, args)
            disp("Check constructor arguments for DIF_LinearDynamics");
            num_vars = this.num_variables;
            assert(isequal(size(args.process_noise_covmat), [num_vars, num_vars]), ...
                "process_noise_covmat is NOT correct size matrix");
        end

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
    end

    methods (Access = protected)
        function predictStateVectorAndCovariance(this)
            Ad = this.discrete_system_matrix;
            Q = this.process_noise_covmat;
            this.state_vector = Ad*this.state_vector;
            this.state_covmat = Ad*this.state_covmat*Ad.' + Q;
        end
    end
end