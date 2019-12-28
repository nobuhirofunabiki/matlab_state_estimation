classdef InformationFilterBase < handle
    properties (SetAccess = immutable)
        num_variables           % Number of variables
    end
    properties (SetAccess = protected)
        info_vector             % xi(k): Information vector
        info_matrix             % Sigma(k): Information matrix
    end
    properties (Abstract = true, SetAccess = protected)
        state_vector            % x(k): State vector
        state_covmat            % P(k): Covariance matrix of state
    end

    methods (Access = protected)
        function obj = InformationFilterBase(args)
            num_vars                    = args.num_variables;
            obj.num_variables           = num_vars;
            obj.info_vector             = zeros(num_vars, 1);
            obj.info_matrix             = zeros(num_vars, num_vars);
        end
    end

    methods (Abstract = true, Access = protected)
        predictStateVectorAndCovariance(this);
    end

    methods (Access = public)
        function executeInformationFilter(this, args)
            this.predictStateVectorAndCovariance();
            this.convertMomentsToInformationForm();
            this.addObservationInformation(...
                args.obs_matrix, args.obs_covmat, args.measures);
            this.convertInformationToMomentsForm();
        end
    end

    methods (Access = protected)
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
    end

    methods (Access = public)
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