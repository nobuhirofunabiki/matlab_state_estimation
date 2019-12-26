classdef UnscentedKalmanFilterBase < handle
    properties (SetAccess = protected)
        alpha
        beta
        gamma
        lambda
        kappa
        state_dim       % state dimunsion (L)
        weights_mean    % weights for means
        weights_cov     % weights for covariances
        sigma_points    % [row: estimated state, column: number of sigma point]
        state_est       % representative estimated state
        state_cov       % state covariance
        X_diff
        Pzz
        Pxz
        Z_diff
        z_pred
    end

    methods
        function obj = UnscentedKalmanFilterBase(args)
            obj.alpha       = args.alpha;
            obj.beta        = args.beta;
            obj.kappa       = args.kappa;
            obj.state_dim   = args.state_dim;
            obj.state_est   = args.state_est;
            obj.state_cov   = args.state_cov;
            obj.lambda      = obj.alpha^2 * (obj.state_dim + obj.kappa) - obj.state_dim;
            gamma_2         = obj.state_dim + obj.lambda;
            obj.gamma       = sqrt(gamma_2);
            % Weights for means
            obj.weights_mean = zeros(1,1+2*obj.state_dim);
            obj.weights_mean(1,1) = obj.lambda/gamma_2;
            for iPoints = 2:1+2*obj.state_dim
                obj.weights_mean(1,iPoints) = 1.0/(2.0*gamma_2);
            end
            % Weights for covariances
            obj.weights_cov = zeros(1,1+2*obj.state_dim);
            obj.weights_cov(1,1) = obj.weights_mean(1,1) ...
                + (1.0 - obj.alpha^2 + obj.beta);
            for iPoints = 2:1+2*obj.state_dim
                obj.weights_cov(1,iPoints) = 1.0/(2.0*gamma_2);
            end
            % Set sigma points
            obj.generateSigmaPoints();
        end

        function generateSigmaPoints(this)
            A = this.gamma * chol(this.state_cov).';
            x = this.state_est;
            Y = x(:,ones(1,numel(x)));
            this.sigma_points = [x Y+A Y-A];
        end

        function executeUnscentedKalmanFilter(this, measurements)
            this.generateSigmaPoints();
            this.predictStates();
            this.generateSigmaPoints();
            this.updateByMeasurement();
            this.Pxz = this.X_diff*diag(this.weights_cov)*this.Z_diff';
            K = this.Pxz*inv(this.Pzz);
            this.state_est = this.state_est + K*(measurements - this.z_pred);
            this.state_cov = this.state_cov - K*this.Pzz*K';
        end

        % Getters ---------------------------------------------
        function output = getStateEstimate(this)
            output = this.state_est;
        end
    end

end