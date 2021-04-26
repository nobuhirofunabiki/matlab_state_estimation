classdef UKF_LinearDynamics < UnscentedKalmanFilterBase
    properties (Abstract = true, SetAccess = protected)
        discrete_system_matrix
        process_noise_covmat        % covariance matrix of process noise
    end
    methods (Access = protected)
        function obj = UKF_LinearDynamics(args)
            obj@UnscentedKalmanFilterBase(args);
        end
    end

    methods (Access = protected)
        function predictStates(this)
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
    end
end