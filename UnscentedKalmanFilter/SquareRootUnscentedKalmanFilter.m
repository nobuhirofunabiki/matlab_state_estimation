classdef SquareRootUnscentedKalmanFilter < UnscentedKalmanFilterBase
    properties (Abstract = true, SetAccess = protected)
        discrete_system_matrix
        process_noise_covmat
        S
        Sz
        sqrtQ
        sqrtR
    end

    methods (Access = public)
        function obj = SquareRootUnscentedKalmanFilter(args)
            obj@UnscentedKalmanFilterBase(args);
        end
    end

    methods (Access = public)
        function executeUnscentedKalmanFilter(this, measurements)
            this.sigma_points = this.generateSigmaPoints();
            this.predictStates();
            this.processMeasurements(measurements);
            this.updateByMeasurements();
        end
    end

    methods (Access = protected)
        function obj = generateSigmaPoints(this)
            A = this.gamma * this.S;
            x = this.state_vector;
            Y = x(:,ones(1,numel(x)));
            obj = [x Y+A Y-A];
        end

        function updateByMeasurements(this)
            this.Pxz = this.X_diff*diag(this.weights_cov)*this.Z_diff';
            this.Pzz = (this.Sz).'*this.Sz;
            K = this.Pxz*inv(this.Pzz);
            this.K = K;
            this.state_vector = this.state_vector + K*(this.measurements - this.z_pred);
            U = K*(this.Sz).';
            for i = 1:size(U,2)
                this.S = cholupdate(this.S, U(:,i), '-');
            end
        end

        function predictStates(this)
            num_sigma_points = size(this.sigma_points, 2);
            state_vector_prior = zeros(size(this.state_vector));
            sigma_point_states = zeros(size(this.sigma_points));
            for iPoints = 1:num_sigma_points
                sigma_point_states(:,iPoints) = ...
                    this.discrete_system_matrix*this.sigma_points(:,iPoints);
                state_vector_prior = state_vector_prior ...
                    + this.weights_mean(1,iPoints)*sigma_point_states(:,iPoints);
            end
            this.state_vector = state_vector_prior;
            this.X_diff = sigma_point_states - state_vector_prior(:,ones(1,num_sigma_points));
            this.S = this.qr_decomposition(this.X_diff, this.sqrtQ);
        end
    end

    methods (Access = protected)
        function obj = qr_decomposition(this, x_diff, square_root_covariance)
            num_sigma_points = size(x_diff,2);
            num_variables = size(x_diff,1);
            x_diff_weighted = zeros(num_variables, num_sigma_points-1);
            for iPoints = 2:num_sigma_points
                x_diff_weighted(:,iPoints-1) = ...
                    sqrt(this.weights_cov(1,iPoints))*x_diff(:,iPoints);
            end
            M = horzcat(x_diff_weighted, square_root_covariance);
            [foo, S] = qr(M.');
            % S = qr(M.');
            S_test = S(1:num_variables,1:num_variables);
            if sign(this.weights_cov(1,1)) == 1
                operation = '+';
            elseif sign(this.weights_cov(1,1)) == -1
                operation = '-';
            else
                error('QR decomposition error');
            end
            obj = cholupdate(S_test, sqrt(norm(this.weights_cov(1,1)))*x_diff(:,1), operation);
        end
    end

end