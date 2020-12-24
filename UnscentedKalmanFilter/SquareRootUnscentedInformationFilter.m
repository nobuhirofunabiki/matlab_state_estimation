classdef SquareRootUnscentedInformationFilter < UnscentedKalmanFilterBase
    properties (Abstract = true, SetAccess = protected)
        discrete_system_matrix
        process_noise_covmat
        info_vector
        sqrt_cov_matrix
        sqrt_info_matrix
        sqrtQ
        sqrtR
    end

    properties (SetAccess = protected)
        functions_
    end

    methods (Access = public)
        function obj = SquareRootUnscentedInformationFilter(args)
            obj@UnscentedKalmanFilterBase(args);
            obj.functions_ = CommonEstimatorFunctions();
        end
    end

    methods (Access = protected)
        function obj = generateSigmaPoints(this)
            sqrt_cov_matrix = transpose(inv(this.sqrt_info_matrix));
            A = this.gamma * sqrt_cov_matrix.';
            x = this.state_vector;
            Y = x(:,ones(1,numel(x)));
            obj = [x Y+A Y-A];
        end

        
        function updateByMeasurements(this)
            this.X_diff = this.sigma_points - this.state_vector;
            this.Pxz = this.X_diff*diag(this.weights_cov)*this.Z_diff';
            Pxz = this.Pxz;
            S_info = this.sqrt_info_matrix;
            Si = inv(this.sqrtR) * transpose(this.Pxz) * transpose(S_info) * S_info;
            HRH = transpose(Si) * Si;
            R = transpose(this.sqrtR) * this.sqrtR;
            % obs_info_vector = (transpose(Si) * inv(this.sqrtR) * (this.measurements - this.z_pred)) + (transpose(Si) * Si * this.state_vector);
            % test1 = Si * inv(this.sqrtR) * (this.measurements - this.z_pred);
            obs_info_vector = transpose(Si)*inv(this.sqrtR)*(this.measurements-this.z_pred) + HRH*this.state_vector;
            Si_square = this.functions_.qr_decomposition_pure(Si.');

            this.info_vector = this.info_vector + obs_info_vector;
            U = Si_square;
            for i = 1:size(U,2)
                this.sqrt_info_matrix = cholupdate(this.sqrt_info_matrix, transpose(U(i,:)), '+');
            end

            this.sqrt_cov_matrix = transpose(inv(this.sqrt_info_matrix));
            this.state_vector = transpose(this.sqrt_cov_matrix) * this.sqrt_cov_matrix * this.info_vector;

            info_vector = this.info_vector;
            info_matrix = transpose(this.sqrt_info_matrix) * this.sqrt_info_matrix;
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
            this.sqrt_cov_matrix = this.functions_.qr_decomposition(this.X_diff, this.sqrtQ, this.weights_cov);
            this.sqrt_info_matrix = transpose(inv(this.sqrt_cov_matrix));
            this.info_vector = transpose(this.sqrt_info_matrix) * this.sqrt_info_matrix * this.state_vector;
            % state_covmat = transpose(this.sqrt_cov_matrix) * this.sqrt_cov_matrix;
            % state_vector = this.state_vector;
        end
    end
end