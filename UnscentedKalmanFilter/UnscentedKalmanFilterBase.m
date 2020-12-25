classdef UnscentedKalmanFilterBase < handle
    properties (SetAccess = protected)
        weights_mean    % weights for means
        weights_cov     % weights for covariances
        sigma_points    % [row: estimated state, column: number of sigma point]
        X_diff
        Pzz
        Pxz
        Z_diff
        z_pred
        measurements
        K
    end
    properties (SetAccess = immutable)
        alpha
        beta
        gamma
        lambda
        kappa
        num_variables   % state dimunsion (L) 
    end
    properties (Abstract = true, SetAccess = protected)
        state_vector    % representative estimated state
        state_covmat    % state covariance
    end

    methods (Access = protected)
        function obj = UnscentedKalmanFilterBase(args)
            obj.alpha           = args.alpha;
            obj.beta            = args.beta;
            obj.kappa           = args.kappa;
            obj.num_variables   = args.num_variables;
            obj.lambda          = obj.alpha^2 * (obj.num_variables + obj.kappa) - obj.num_variables;
            gamma_2             = obj.num_variables + obj.lambda;
            obj.gamma           = sqrt(gamma_2);
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

    methods (Abstract, Access = protected)
        predictStates(this);
        processMeasurements(this);
    end

    methods (Access = public)
        function executeUnscentedKalmanFilter(this, measurements)
            this.sigma_points = this.generateSigmaPoints();
            this.predictStates();
            this.sigma_points = this.generateSigmaPoints();
            this.processMeasurements(measurements);
            this.updateByMeasurements();
        end

        % Getters ---------------------------------------------
        function output = getStateVector(this)
            output = this.state_vector;
        end

        function output = getStateCovarianceMatrix(this)
            output = this.state_covmat;
        end
    end

    methods (Access = protected)
        function obj = generateSigmaPoints(this)
            if (isdiag(this.state_covmat) == false)
                this.state_covmat = 0.5*(this.state_covmat + (this.state_covmat)');
            end
            A = this.gamma * chol(this.state_covmat, 'lower');
            x = this.state_vector;
            Y = x(:,ones(1,numel(x)));
            obj = [x Y+A Y-A];
        end

        function updateByMeasurements(this)
            this.X_diff = this.sigma_points - this.state_vector;
            this.Pxz = this.X_diff*diag(this.weights_cov)*this.Z_diff';
            K = this.Pxz*inv(this.Pzz);
            this.K = K;
            this.state_vector = this.state_vector + K*(this.measurements - this.z_pred);
            this.state_covmat = this.state_covmat - K*this.Pzz*K';
        end
    end

end