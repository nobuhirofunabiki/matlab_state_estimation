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
            temp = transpose(inv(this.sqrt_info_matrix));
            sqrt_cov_matrix = this.functions_.qr_decomposition_pure(temp);
            A = this.gamma * sqrt_cov_matrix;
            x = this.state_vector;
            Y = x(:,ones(1,numel(x)));
            obj = [x Y+A Y-A];
        end

        
        function updateByMeasurements(this)
            this.X_diff = this.sigma_points - this.state_vector;
            this.Pxz = this.X_diff*diag(this.weights_cov)*this.Z_diff'; % PxzはSRUKFと等しい
            S_info = this.sqrt_info_matrix; % info_matrix (=S_info*S_info.') はSRUKFと等しい
            Si = S_info*transpose(S_info)*this.Pxz/this.sqrtR;

            % test1とtest2の結果が同じになる...何かが根本的に間違っていそう
            % ただしtest2をSRUKF側で計算した結果とは等しくなっている...
            test1 = Si*transpose(Si);
            P = this.sqrt_cov_matrix * transpose(this.sqrt_cov_matrix);
            Pxz_2 = this.Pxz;
            H = transpose(Pxz_2)*inv(P);
            R = this.sqrtR * transpose(this.sqrtR);
            test2 = H.'/R*H;

            % this.measurements, this.z_pred, this.state_vectorいずれもSRUKFと等しい
            % obs_info_vector1の値がEIFと大きく異なる
            % EIFと大きく異なったのはthis.z_predの値、obs_info_vector1_confは全く同じではないがEIFとほぼ等しかった...
            obs_info_vector1_cof = Si * transpose(inv(this.sqrtR));
            obs_info_vector1 = Si * transpose(inv(this.sqrtR)) * (this.measurements - this.z_pred);
            obs_info_vector2 = Si * transpose(Si) * this.state_vector;
            obs_info_vector = obs_info_vector1 + obs_info_vector2;
            % obs_info_vector = Si * transpose(inv(this.sqrtR)) * (this.measurements - this.z_pred) + Si * transpose(Si) * this.state_vector;
            this.info_vector = this.info_vector + obs_info_vector;
            
            U = this.functions_.qr_decomposition_pure(Si);
            hoge = transpose(this.sqrt_info_matrix);
            for i = 1:size(U,2)
                % this.sqrt_info_matrix = cholupdate(this.sqrt_info_matrix, transpose(U(i,:)), '+');
                hoge = cholupdate(hoge, transpose(U(i,:)), '+');
            end
            this.sqrt_info_matrix = hoge.';

            temp = transpose(inv(this.sqrt_info_matrix));
            this.sqrt_cov_matrix= this.functions_.qr_decomposition_pure(temp);
            this.state_vector = this.sqrt_cov_matrix * transpose(this.sqrt_cov_matrix) * this.info_vector;
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
            hoge = transpose(inv(this.sqrt_cov_matrix));
            this.sqrt_info_matrix = this.functions_.qr_decomposition_pure(hoge);
            this.info_vector = this.sqrt_info_matrix * transpose(this.sqrt_info_matrix) * this.state_vector;
        end
    end
end