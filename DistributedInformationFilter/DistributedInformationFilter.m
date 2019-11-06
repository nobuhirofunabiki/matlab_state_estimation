classdef DistributedInformationFilter < handle
    properties (SetAccess = protected)
        agent_id
        state_vector                % x(k)
        state_covmat                % P(k)
        info_vector                 % z(k)
        info_matrix                 % Z(k)
        number_agents               % N
        obs_info_vector             % i(k)
        obs_info_vector_prev        % i(k-1)
        obs_info_vector_joint       % u(k)
        obs_info_vector_joint_prev  % u(k-1)
        obs_info_matrix             % I(k)
        obs_info_matrix_prev        % I(k-1)
        obs_info_matrix_joint       % U(k)
        obs_info_matrix_joint_prev  % U(k-1)
        out_info_vector_pool        % 
        out_info_matrix_pool        %
        system_matrix               % Ad
    end
    methods
        function obj = DistributedInformationFilter(args)
            NUM_VAR = args.number_variables;
            obj.agent_id                    = args.agent_id;
            obj.number_agents               = args.number_agents;
            obj.state_vector                = zeros(NUM_VAR, 1);
            obj.state_covmat                = zeros(NUM_VAR, NUM_VAR);
            obj.info_vector                 = zeros(NUM_VAR, 1);
            obj.info_matrix                 = zeros(NUM_VAR, NUM_VAR);
            obj.obs_info_vector             = zeros(NUM_VAR, 1);
            obj.obs_info_vector_prev        = zeros(NUM_VAR, 1);
            obj.obs_info_vector_joint       = zeros(NUM_VAR, 1);
            obj.obs_info_vector_joint_prev  = zeros(NUM_VAR, 1);
            obj.obs_info_matrix             = zeros(NUM_VAR, NUM_VAR);
            obj.obs_info_matrix_prev        = zeros(NUM_VAR, NUM_VAR);
            obj.obs_info_matrix_joint       = zeros(NUM_VAR, NUM_VAR);
            obj.obs_info_matrix_joint_prev  = zeros(NUM_VAR, NUM_VAR);
            obj.out_info_vector_pool        = zeros(NUM_VAR, 1);
            obj.out_info_matrix_pool        = zeros(NUM_VAR, NUM_VAR);
            obj.system_matrix               = zeros(NUM_VAR, NUM_VAR);
        end
        function propagateCovarianceMatrix(this)
            Ad = this.system_matrix;
            % TODO: Add control and disturbance noises
            this.state_covmat = Ad*this.state_covmat*Ad.';
        end
        function updateEstimatorStatus(this)
            this.obs_info_vector_prev = this.obs_info_vector;
            this.obs_info_vector_joint_prev = this.obs_info_vector_joint;
            this.obs_info_matrix_prev = this.obs_info_matrix;
            this.obs_info_matrix_joint_prev = this.obs_info_matrix_joint;
        end
        function integrateMultiSourceInformation(this)
            this.obs_info_vector_joint = ...
                this.obs_info_vector - this.obs_info_vector_prev ...
                + this.out_info_vector_pool;
            this.obs_info_matrix_joint = ...
                this.obs_info_matrix - this.obs_info_matrix_prev ...
                + this.out_info_matrix_pool;
        end
        function computePosteriorPdf(this)
            N = this.number_agents;
            this.info_matrix = inv(this.state_covmat);
            this.info_vector = this.info_matrix * this.state_vector;
            this.info_vector = this.info_vector + N * this.obs_info_vector_joint;
            this.info_matrix = this.info_matrix + N * this.obs_info_matrix_joint;
            this.state_covmat = inv(this.info_matrix);
            this.state_vector = this.state_covmat * this.info_vector;
        end

        function addObservationInformation(this, obs_matrix, obs_covmat, obs_measured)
            H = obs_matrix;
            R = obs_covmat;
            y = obs_measured;
            this.obs_info_vector = this.obs_info_vector + H.'/R*y;
            this.obs_info_matrix = this.obs_info_matrix + H.'/R*H;
        end
        function resetObservationInformation(this)
            this.obs_info_vector = zeros(size(this.obs_info_vector));
            this.obs_info_matrix = zeros(size(this.obs_info_matrix));
        end
        function clearOutSourceInformationPool(this)
            this.out_info_vector_pool = zeros(size(this.out_info_vector_pool));
            this.out_info_matrix_pool = zeros(size(this.out_info_matrix_pool));
        end
        
        % Setters ---------------------------------------------------

        function addOutSourceInformationIntoPool(this, info_vector, info_matrix)
            this.out_info_vector_pool = this.out_info_vector_pool + info_vector;
            this.out_info_matrix_pool = this.out_info_matrix_pool + info_matrix;
        end
        function setStateCovarianceMatrix(this, args)
            P = zeros(size(this.state_covmat));
            % TODO
            num_dims = 2;
            for iAgents = 1:this.number_agents
                for iDims = 1:num_dims
                    P(2*num_dims*(iAgents-1)+iDims, 2*num_dims*(iAgents-1)+iDims) = args.position_sigma^2;
                    P(2*num_dims*(iAgents-1)+num_dims+iDims, 2*num_dims*(iAgents-1)+num_dims+iDims) = args.velocity_sigma^2;
                end
            end
            this.state_covmat = P;
        end
        function setDiscreteSystemMatrix(this, discrete_system_matrix)
            NUM_AGENTS = this.number_agents;
            % TODO
            num_dims = 2;
            num_vars = 2*num_dims;
            A = zeros(num_vars*NUM_AGENTS, num_vars*NUM_AGENTS);
            for iAgents = 1:NUM_AGENTS
                A(1+num_vars*(iAgents-1):num_vars*iAgents,...
                    1+num_vars*(iAgents-1):num_vars*iAgents) = discrete_system_matrix;
            end
            this.system_matrix = A;
        end
        function setEstimatedVariable(this, index_begin, index_end, arg_variable)
            this.state_vector(index_begin:index_end, 1) = arg_variable;
        end

        % Getters -------------------------------------------------------

        function output = getStateVector(this)
            output = this.state_vector;
        end
        function output = getStateCovarianceMatrix(this)
            output = this.state_covmat;
        end
        function output = getCurrentJointInformationVector(this)
            output = this.obs_info_vector_joint;
        end
        function output = getPreviousJointInformationVector(this)
            output = this.obs_info_vector_joint_prev;
        end
        function output = getCurrentJointInformationMatrix(this)
            output = this.obs_info_matrix_joint;
        end
        function output = getPreviousJointInformationMatrix(this)
            output = this.obs_info_matrix_joint_prev;
        end
    end
end