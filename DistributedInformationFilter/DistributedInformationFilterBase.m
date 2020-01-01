classdef DistributedInformationFilterBase < handle
    properties (SetAccess = protected)
        agent_id
        num_variables
        state_vector                % x(k)
        state_covmat                % P(k)
        info_vector                 % z(k)
        info_matrix                 % Z(k)
        num_agents                  % N
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
    end

    methods (Access = protected)
        function obj = DistributedInformationFilterBase(args)
            obj.checkConstructorArguments(args);
            NUM_VAR = args.number_variables;
            obj.num_variables               = NUM_VAR;
            obj.agent_id                    = args.agent_id;
            obj.num_agents                  = args.num_agents;
            obj.state_vector                = args.state_vector;
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
        end
    end

    methods (Access = private)
        function checkConstructorArguments(this, args)
            disp("Check constructor arguments for DistributedInformationFilterBase");
            assert(isequal(size(args.state_vector), [args.number_variables, 1]), ...
                "state_vector is NOT correct size vector");
        end
    end

    methods (Abstract, Access = protected)
        processMeasurements(this);
        predictStateVectorAndCovariance(this);
    end

    methods (Access = public)
        function executeFiltering(this, measures, adjacent_matrix, ...
            outsource_info_vector, outsource_info_matrix)
            this.predictStateVectorAndCovariance();
            this.resetObservationInformation();
            this.processMeasurements(measures, adjacent_matrix);
            this.clearOutSourceInformationPool();
            this.addOutSourceInformationIntoPool(...
                outsource_info_vector, outsource_info_matrix);
            this.integrateMultiSourceInformation();
            this.computePosteriorPdf();
        end

        function updateEstimatorStatus(this)
            this.obs_info_vector_prev = this.obs_info_vector;
            this.obs_info_vector_joint_prev = this.obs_info_vector_joint;
            this.obs_info_matrix_prev = this.obs_info_matrix;
            this.obs_info_matrix_joint_prev = this.obs_info_matrix_joint;
        end
    end

    methods (Access = protected)
        function resetObservationInformation(this)
            this.obs_info_vector = zeros(size(this.obs_info_vector));
            this.obs_info_matrix = zeros(size(this.obs_info_matrix));
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
            N = this.num_agents;
            this.info_matrix = inv(this.state_covmat);
            this.info_vector = this.info_matrix * this.state_vector;
            this.info_vector = this.info_vector + N * this.obs_info_vector_joint;
            this.info_matrix = this.info_matrix + N * this.obs_info_matrix_joint;
            this.state_covmat = inv(this.info_matrix);
            this.state_vector = this.state_covmat * this.info_vector;
        end

        function clearOutSourceInformationPool(this)
            this.out_info_vector_pool = zeros(size(this.out_info_vector_pool));
            this.out_info_matrix_pool = zeros(size(this.out_info_matrix_pool));
        end

        function addOutSourceInformationIntoPool(this, info_vector, info_matrix)
            this.out_info_vector_pool = this.out_info_vector_pool + info_vector;
            this.out_info_matrix_pool = this.out_info_matrix_pool + info_matrix;
        end
    end

    methods (Access = protected)
        function addObservationInformation(this, obs_matrix, obs_covmat, obs_measured)
            H = obs_matrix;
            R = obs_covmat;
            y = obs_measured;
            this.obs_info_vector = this.obs_info_vector + H.'/R*y;
            this.obs_info_matrix = this.obs_info_matrix + H.'/R*H;
        end
        
        % Setters ---------------------------------------------------
        
        function setEstimatedVariable(this, index_begin, index_end, arg_variable)
            this.state_vector(index_begin:index_end, 1) = arg_variable;
        end
    end

    methods (Access = public)
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