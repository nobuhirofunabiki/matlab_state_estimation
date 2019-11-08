classdef InformationFilter < handle
    properties (SetAccess = protected)
        num_variables       % Number of variables
        num_dims            % Number of dimensions (2D or 3D)
        num_agents          % Number of agents
        state_vector        % x(k): State vector
        state_covmat        % P(k): Covariance matrix of state
        info_vector         % xi(k): Information vector
        info_matrix         % Sigma(k): Information matrix
        system_matrix       % Ad(k): System matrix
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
            obj.system_matrix   = zeros(num_vars, num_vars);
        end
        function propagateCovarianceMatrix(this)
            Ad = this.system_matrix;
            % TODO: Add control and disturbance noises
            this.state_covmat = Ad*this.state_covmat*Ad.';
        end
        function addObservationInformation(this, obs_matrix, obs_covmat, obs_measured)
            H = obs_matrix;
            R = obs_covmat;
            y = obs_measured;
            this.info_vector = this.info_vector + H.'/R*y;
            this.info_matrix = this.info_matrix + H.'/R*H;
        end
        function convertMomentsToInformationForm(this)
            this.info_matrix = inv(this.state_covmat);
            this.info_vector =  this.info_matrix * this.state_vector;
        end

        % Setters
        function setEstimatedVariable(this, index_begin, index_end, arg_variable)
            this.state_vector(index_begin:index_end, 1) = arg_variable;
        end
    end
end