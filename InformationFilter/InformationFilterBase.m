classdef InformationFilterBase < handle
    properties (SetAccess = protected)
        num_variables           % Number of variables
        num_dimensions          % Number of dimensions (2D or 3D)
        state_vector            % x(k): State vector
        state_covmat            % P(k): Covariance matrix of state
        process_noise_covmat    % Q(k): Covariance matrix of process noise
        info_vector             % xi(k): Information vector
        info_matrix             % Sigma(k): Information matrix
        discrete_system_matrix  % Ad(k): Discrete system matrix
    end

    methods
        function obj = InformationFilterBase(args)
            num_vars                    = args.num_variables;
            obj.num_variables           = num_vars;
            obj.num_dimensions          = args.num_dimensions;
            obj.info_vector             = zeros(num_vars, 1);
            obj.info_matrix             = zeros(num_vars, num_vars);
        end

        function executeInformationFilter(this, args)
            this.predictStateVectorAndCovariance();
            this.convertMomentsToInformationForm();
            this.addObservationInformation(...
                args.obs_matrix, args.obs_covmat, args.measures);
            this.convertInformationToMomentsForm();
        end

        function predictStateVectorAndCovariance(this)
            Ad = this.discrete_system_matrix;
            Q = this.process_noise_covmat;
            this.state_vector = Ad*this.state_vector;
            this.state_covmat = Ad*this.state_covmat*Ad.' + Q;
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

        function setStateCovarianceMatrixAll(this, state_covmat)
            this.state_covmat = state_covmat;
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