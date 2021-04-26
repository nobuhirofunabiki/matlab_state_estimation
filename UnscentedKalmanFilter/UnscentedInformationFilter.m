classdef UnscentedInformationFilter < UnscentedKalmanFilterBase
    properties (Abstract = true, SetAccess = protected)
        discrete_system_matrix
        process_noise_covmat
        info_vector
        info_matrix
        R
    end

    methods (Access = public)
        function obj = UnscentedInformationFilter(args)
            obj@UnscentedKalmanFilterBase(args);
        end
    end

    methods (Access = protected)
        function updateByMeasurements(this)
            this.X_diff = this.sigma_points - this.state_vector;
            this.Pxz = this.X_diff*diag(this.weights_cov)*this.Z_diff';
            R = this.R;
            Pxz = this.Pxz;
            Y = this.info_matrix;
            y = this.info_vector;
            info_contribution_vector = Y*Pxz/R*(this.measurements-this.z_pred+Pxz.'*y);
            info_contribution_matrix = Y*Pxz/R*Pxz.'*Y.';
            this.info_vector = this.info_vector + info_contribution_vector;
            this.info_matrix = this.info_matrix + info_contribution_matrix;
            this.state_covmat = inv(this.info_matrix);
            this.state_vector = this.state_covmat * this.info_vector;
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
            this.state_covmat = this.X_diff*diag(this.weights_cov)*(this.X_diff).' + this.process_noise_covmat;
            this.info_matrix = inv(this.state_covmat);
            this.info_vector = this.info_matrix * this.state_vector;
        end
    end
end