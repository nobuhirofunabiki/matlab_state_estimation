classdef UKF_LinearDynamics < UnscentedKalmanFilterBase
    properties (SetAccess = protected)
        discrete_system_matrix
        process_noise_covmat 
    end

    methods (Access = protected)
        function obj = UKF_LinearDynamics(args)
            obj@UnscentedKalmanFilterBase(args);
            obj.discrete_system_matrix = args.discrete_system_matrix;
            obj.process_noise_covmat = args.process_noise_covmat;
        end
    end

    methods (Access = protected)
        function predictStates(this)
            num_sigma_points = size(this.sigma_points, 2);
            state_est_prior = zeros(size(this.state_est));
            sigma_point_states = zeros(size(this.sigma_points));
            for iPoints = 1:num_sigma_points
                sigma_point_states(:,iPoints) = this.discrete_system_matrix*this.sigma_points(:,iPoints);
                state_est_prior = state_est_prior ...
                    + this.weights_mean(1,iPoints)*sigma_point_states(:,iPoints);
            end
            this.state_est = state_est_prior;
            this.X_diff = sigma_point_states - state_est_prior(:,ones(1,num_sigma_points));
            this.state_cov = this.X_diff*diag(this.weights_cov)*(this.X_diff)' + this.process_noise_covmat;
        end
    end
end