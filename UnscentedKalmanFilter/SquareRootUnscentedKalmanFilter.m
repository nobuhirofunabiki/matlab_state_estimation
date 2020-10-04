classdef SquareRootUnscentedKalmanFilter < UnscentedKalmanFilterBase
    properties (Abstract = true, SetAccess = protected)
        discrete_system_matrix
        process_noise_covmat
        S
        Sz
        sqrtQ
        sqrtR
    end

    properties (SetAccess = protected)
        functions_
    end

    methods (Access = public)
        function obj = SquareRootUnscentedKalmanFilter(args)
            obj@UnscentedKalmanFilterBase(args);
            obj.functions_ = CommonEstimatorFunctions();
        end
    end

    methods (Access = protected)
        function obj = generateSigmaPoints(this)
            A = this.gamma * (this.S).';
            x = this.state_vector;
            Y = x(:,ones(1,numel(x)));
            obj = [x Y+A Y-A];
        end

        function updateByMeasurements(this)
            this.X_diff = this.sigma_points - this.state_vector;
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
            this.S = this.functions_.qr_decomposition(this.X_diff, this.sqrtQ, this.weights_cov);
        end
    end

end