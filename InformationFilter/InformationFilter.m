classdef InformationFilter < handle
    properties (SetAccess = protected)
        num_variables       % Number of variables
        num_dims            % Number of dimensions (2D or 3D)
        num_agents          % Number of agents
        state_vector        % x(k): State vector
        state_covmat        % P(k): Covariance matrix of state
        info_vector         % xi(k): Information vector
        info_matrix         % Sigma(k): Information matrix
        discrete_system_matrix       % Ad(k): System matrix
    end
    methods
        function obj = InformationFilter(args)
            num_vars            = args.num_variables;
            obj.num_variables   = num_vars;
            obj.num_dims        = args.num_dims;
            obj.num_agents      = args.num_agents;
            obj.state_vector    = zeros(num_vars, 1);
            obj.state_covmat    = zeros(num_vars, num_vars);
            obj.info_vector     = zeros(num_vars, 1);
            obj.info_matrix     = zeros(num_vars, num_vars);
            obj.discrete_system_matrix   = zeros(num_vars, num_vars);
        end
        function executeInformationFilter(this, args)
            this.setStateVectorAll(args.prior_state);
            this.convertMomentsToInformationForm();
            this.addObservationInformation(...
                args.obs_matrix, args.obs_covmat, args.measures);
            this.convertInformationToMomentsForm();
        end
        function propagateCovarianceMatrix(this)
            Ad = this.discrete_system_matrix;
            % TODO: Add control and disturbance noises
            this.state_covmat = Ad*this.state_covmat*Ad.';
        end
        function addObservationInformation(this, obs_matrix, obs_covmat, measures)
            H = obs_matrix;
            R = obs_covmat;
            y = measures;
            this.info_vector = this.info_vector + H.'/R*y;
            this.info_matrix = this.info_matrix + H.'/R*H;
        end
        function convertMomentsToInformationForm(this)
            this.info_matrix = inv(this.state_covmat);
            this.info_vector =  this.info_matrix * this.state_vector;
        end
        function convertInformationToMomentsForm(this)
            this.state_covmat = inv(this.info_matrix);
            this.state_vector = this.state_covmat * this.info_vector;
        end

        % Setters-------------------------------------------------------
        function setEstimatedVariable(this, index_begin, index_end, arg_variable)
            this.state_vector(index_begin:index_end, 1) = arg_variable;
        end
        function setStateVectorAll(this, arg)
            this.state_vector = arg;
        end
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
        function setStateCovarianceMatrixAll(this, state_covmat)
            this.state_covmat = state_covmat;
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

        % Getters -------------------------------------------------------
        function output = getStateVector(this)
            output = this.state_vector;
        end
        function output = getStateCovarianceMatrix(this)
            output = this.state_covmat;
        end
    end
end