classdef DecentralizedUnscentedInformationFilterBase < DecentralizedInformationFilterBase
    properties (Abstract = true, SetAccess = immutable)
        num_dimensions
        discrete_system_matrix
        process_noise_covmat
    end
    properties (Abstract = true, SetAccess = protected)
        state_vector
        state_covmat
    end
    properties (SetAccess = immutable)
        alpha
        beta
        gamma
        lambda
        kappa
        weights_mean    % weights for means
        weights_cov     % weights for covariances
    end
    properties (SetAccess = protected)
        sigma_points    % [row: estimated state, column: number of sigma point]
        X_diff
    end

    methods (Access = protected)
        function obj = DecentralizedUnscentedInformationFilterBase(args)
            obj@DecentralizedInformationFilterBase(args);
            obj.process_noise_covmat    = args.process_noise_covmat;
            obj.alpha                   = args.alpha;
            obj.beta                    = args.beta;
            obj.kappa                   = args.kappa;
            obj.lambda                  = obj.alpha^2 * (obj.num_variables + obj.kappa) - obj.num_variables;
            gamma_2                     = obj.num_variables + obj.lambda;
            obj.gamma                   = sqrt(gamma_2);
            % Weights for means
            obj.weights_mean = zeros(1,1+2*obj.num_variables);
            obj.weights_mean(1,1) = obj.lambda/gamma_2;
            for iPoints = 2:1+2*obj.num_variables
                obj.weights_mean(1,iPoints) = 1.0/(2.0*gamma_2);
            end
            % Weights for covariances
            obj.weights_cov = zeros(1,1+2*obj.num_variables);
            obj.weights_cov(1,1) = obj.weights_mean(1,1) ...
                + (1.0 - obj.alpha^2 + obj.beta);
            for iPoints = 2:1+2*obj.num_variables
                obj.weights_cov(1,iPoints) = 1.0/(2.0*gamma_2);
            end
        end
    end

    methods (Access = protected)
        function predictStateVectorAndCovariance(this)
            num_sigma_points = size(this.sigma_points, 2);
            state_vector_prior = zeros(size(this.state_vector));
            sigma_point_states = zeros(size(this.sigma_points));
            for iPoints = 1:num_sigma_points
                sigma_point_states(:,iPoints) = this.discrete_system_matrix*this.sigma_points(:,iPoints);
                state_vector_prior = state_vector_prior ...
                    + this.weights_mean(1,iPoints)*sigma_point_states(:,iPoints);
            end
            this.state_vector = state_vector_prior;
            this.X_diff = sigma_point_states - state_vector_prior(:,ones(1,num_sigma_points));
            this.state_covmat = this.X_diff*diag(this.weights_cov)*(this.X_diff).' + this.process_noise_covmat;
        end

        function output = generateSigmaPoints(this, mean_, covariance_)
            if (isdiag(covariance_) == false)
                covariance_ = 0.5*(covariance_ + (covariance_)');
                this.state_covmat = covariance_;
            end
            Y = mean_(:,ones(1,numel(mean_)));
            S = chol(covariance_);
            A = this.gamma * S.';
            output = [mean_ Y+A Y-A];
        end
    end

end