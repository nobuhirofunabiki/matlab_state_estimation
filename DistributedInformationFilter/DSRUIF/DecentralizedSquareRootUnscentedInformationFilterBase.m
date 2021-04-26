classdef DecentralizedSquareRootUnscentedInformationFilterBase < DecentralizedInformationFilterBase
    properties (Abstract = true, SetAccess = protected)
        sqrtQ
        sqrt_info_matrix % S(k)
    end
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
        sqrtR
        sqrt_obs_info_matrix                % S_i(k)
        sqrt_obs_info_matrix_prev           % S_i(k-1)
        sqrt_obs_info_matrix_joint          % Us(k)
        sqrt_obs_info_matrix_joint_prev     % Us(k-1)
        sqrt_out_info_matrix_pool           %
        step = 0;
    end
    properties (SetAccess = protected)
        functions_ = CommonEstimatorFunctions()
    end

    methods (Access = protected)
        function obj = DecentralizedSquareRootUnscentedInformationFilterBase(args)
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
            % Square-root information filter
            NUM_VAR                                 = args.number_variables;
            obj.sqrt_info_matrix                    = zeros(NUM_VAR, NUM_VAR);
            obj.sqrt_obs_info_matrix                = zeros(NUM_VAR, NUM_VAR);
            obj.sqrt_obs_info_matrix_prev           = zeros(NUM_VAR, NUM_VAR);
            obj.sqrt_obs_info_matrix_joint          = zeros(NUM_VAR, NUM_VAR);
            obj.sqrt_obs_info_matrix_joint_prev     = zeros(NUM_VAR, NUM_VAR);
            obj.sqrt_out_info_matrix_pool           = zeros(NUM_VAR, NUM_VAR);
        end
    end

    methods (Access = public)
        function executeFiltering(this, measures, adjacent_matrix, ...
            outsource_info_vector, outsource_info_matrix)
            this.predictStateVectorAndCovariance();
            this.resetObservationInformation();
            this.processMeasurements(measures, adjacent_matrix);
            this.clearOutSourceInformationPool();
            this.addOutSourceInformationIntoPool(...
                outsource_info_vector, outsource_info_matrix);
            this.integrateMultiSourceInformation();
            this.computePosteriorPdf();
        end

        function updateEstimatorStatus(this)
            this.obs_info_vector_prev = this.obs_info_vector;
            this.obs_info_vector_joint_prev = this.obs_info_vector_joint;
            this.sqrt_obs_info_matrix_prev = this.sqrt_obs_info_matrix;
            this.sqrt_obs_info_matrix_joint_prev = this.sqrt_obs_info_matrix_joint;
        end

        function output = getCurrentJointSquareRootInformationMatrix(this)
            output = this.sqrt_obs_info_matrix_joint;
        end
        function output = getPreviousJointSquareRootInformationMatrix(this)
            output = this.sqrt_obs_info_matrix_joint_prev;
        end
    end

    methods (Access = protected)
        function output = generateSigmaPoints(this, mean_, sqrt_covariance_)
            Y = mean_(:,ones(1,numel(mean_)));
            A = this.gamma * sqrt_covariance_.';
            output = [mean_ Y+A Y-A];
        end

        function predictStateVectorAndCovariance(this)
            sqrt_cov_matrix = transpose(inv(this.sqrt_info_matrix));
            this.sigma_points = this.generateSigmaPoints(this.state_vector, sqrt_cov_matrix);
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
            S = this.functions_.qr_decomposition(this.X_diff, this.sqrtQ, this.weights_cov);
            this.sqrt_info_matrix = transpose(inv(S));
            this.info_vector = transpose(this.sqrt_info_matrix) * this.sqrt_info_matrix * this.state_vector;
            info_matrix = transpose(this.sqrt_info_matrix) * this.sqrt_info_matrix;
        end

        function resetObservationInformation(this)
            this.obs_info_vector = zeros(size(this.obs_info_vector));
            this.sqrt_obs_info_matrix = zeros(size(this.sqrt_obs_info_matrix));
        end

        function clearOutSourceInformationPool(this)
            this.out_info_vector_pool = zeros(size(this.out_info_vector_pool));
            this.sqrt_out_info_matrix_pool = zeros(size(this.sqrt_out_info_matrix_pool));
        end

        function addOutSourceInformationIntoPool(this, info_vector, sqrt_info_matrix)
            this.out_info_vector_pool = this.out_info_vector_pool + info_vector;
            % this.sqrt_out_info_matrix_pool = this.sqrt_out_info_matrix_pool + sqrt_info_matrix;
            hoge = this.sqrt_out_info_matrix_pool;
            for i = 1:size(sqrt_info_matrix, 1)
                % hoge = cholupdate(hoge, transpose(sqrt_info_matrix(i,:)), '+');
                hoge = cholupdate(hoge, sqrt_info_matrix(:,i), '+');
            end
            this.sqrt_out_info_matrix_pool = hoge;
        end

        function integrateMultiSourceInformation(this)
            this.step = this.step + 1;
            if (this.step == 1)
                % this.obs_info_vector_joint = this.obs_info_vector;
                % this.sqrt_obs_info_matrix_joint = this.sqrt_obs_info_matrix;
            else
                % Information contribution vector
                % this.obs_info_vector_joint = ...
                %     this.obs_info_vector - this.obs_info_vector_prev ...
                %     + this.out_info_vector_pool;
                this.obs_info_vector_joint = ...
                    this.obs_info_vector_prev + this.out_info_vector_pool;
                % Square-root information contribution matrix
                % this.sqrt_obs_info_matrix_joint = this.sqrt_obs_info_matrix;
                this.sqrt_obs_info_matrix_joint = zeros(size(this.sqrt_obs_info_matrix_joint));
                hoge = this.sqrt_obs_info_matrix_joint;
                % obs_info_matrix = hoge.'*hoge;
                
                U2 = this.sqrt_out_info_matrix_pool;
                % test5 = U2.'*U2;
                for i = 1:size(U2,1)
                    % hoge = cholupdate(hoge, transpose(U2(i,:)), '+');
                    hoge = cholupdate(hoge, U2(:,i), '+');
                end

                % test3 = hoge.'*hoge + U2.'*U2;

                % test2 = hoge.'*hoge;

                % こいつが上手くいかない...
                % mainから読み込んでいるout_info_poolがduifと一致していることは確認済み
                % sqrt_obs_info_matrix_prevがduifと一致していることは確認済み
                % U1 = this.sqrt_obs_info_matrix_prev;
                % obs_info_matrix_prev = U1.'*U1;
                % for i = 1:size(U1,1)
                %     % hoge = cholupdate(hoge, transpose(U1(i,:)), '-');
                %     hoge = cholupdate(hoge, U1(:,i), '-');
                % end
                % test1 = hoge.'*hoge - U1.'*U1;
                % hoge = chol(test1);

                this.sqrt_obs_info_matrix_joint = hoge;
            end
            obs_info_matrix_joint = transpose(this.sqrt_obs_info_matrix_joint) * this.sqrt_obs_info_matrix_joint; % Delete this later
        end

        function computePosteriorPdf(this)
            N = this.num_agents;
            U = sqrt(N)*this.sqrt_obs_info_matrix_joint;
            hoge = this.sqrt_info_matrix;
            for i = 1:size(U,1)
                % hoge = cholupdate(hoge, transpose(U(i,:)), '+');
                hoge = cholupdate(hoge, U(:,i), '+');
            end

            this.sqrt_info_matrix = hoge;

            info_matrix = (this.sqrt_info_matrix).'*this.sqrt_info_matrix;

            this.info_vector = this.info_vector + N * this.obs_info_vector_joint;

            this.state_covmat = inv(info_matrix);
            this.state_vector = this.state_covmat * this.info_vector;

            this.obs_info_vector_joint = this.obs_info_vector_joint + this.obs_info_vector;
            U2 = this.sqrt_out_info_matrix_pool;
            for i = 1:size(U2,1)
                this.sqrt_obs_info_matrix_joint = cholupdate(this.sqrt_obs_info_matrix_joint, U2(:,i), '+');
            end
        end
    end
end