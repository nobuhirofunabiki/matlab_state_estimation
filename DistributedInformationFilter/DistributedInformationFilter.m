classdef DistributedInformationFilter < handle
    properties (SetAccess = protected)
        state_vector                % x(k)
        state_covmat                % P(k)
        info_vector                 % z(k)
        info_matrix                 % Z(k)
        obs_matrix                  % H(k)
        obs_covmat                  % R(k)
        obs_info_vector             % i(k)
        obs_info_vector_prev        % i(k-1)
        obs_info_vector_joint       % u(k)
        obs_info_vector_joint_prev  % u(k-1)
        obs_info_matrix             % I(k)
        obs_info_matrix_prev        % I(k-1)
        obs_info_matrix_joint       % U(k)
        obs_info_matrix_joint_prev  % U(k-1)
        out_info_vector_pool        % 
        out_info_matrix_pool        %
        obs_mesured                 % y(k)
        Ad
    end
    methods
        function obj = DistributedInformationFilter(args)
            
        end
        function propagateCovarianceMatrix(this)
            Ad = this.Ad;
            % TODO: Add control and disturbance noises
            this.state_covmat = Ad*this.state_covmat*Ad.';
        end
        function calculateObservationInformation(this)
            H = this.obs_matrix;
            R = this.obs_covmat;
            this.obs_info_vector = H.'/R*y;
            this.obs_info_matrix = H.'/R*H;
        end
        function integrateMultiSourceInformation(this)
            this.obs_info_vector_joint = ...
                this.obs_info_vector - this.obs_info_vector_prev ...
                + this.out_info_vector_pool;
            this.obs_info_matrix_joint = ...
                this.obs_info_matrix - this.obs_info_matrix_prev ...
                + this.out_info_matrix_pool;
        end
        function computePosteriorPdf(this)
            this.info_vector = this.info_vector + N * this.obs_info_vector_joint;
            this.info_matrix = this.info_matrix + N * this.obs_info_matrix_joint;
            this.state_covmat = inv(this.info_matrix);
            this.state_vector = this.state_covmat + this.info_vector;

        end
        
        % Setters ---------------------------------------------------

        function addOutSourceInformationIntoPool(this, info_vector, info_matrix)
            this.out_info_vector_pool = this.out_info_vector_pool + info_vector;
            this.out_info_matrix_pool = this.out_info_matrix_pool + info_matrix;
        end
        function setMeasurementData(this, obs_measured)
            this.obs_measured  = obs_measured;
        end
        function setObservationMatrix(this, pos_i, pos_j, agent_id_i, agent_id_j, obs_index)
            dist = norm(pos_i - pos_j);
            DIM  = numel(pos_i);
            for iDim = 1:DIM
                % Agent ID: i
                this.obs_matrix(obs_index, 2*DIM*(agent_id_i-1)+iDim)     = (pos_i(iDim,1)-pos_j(iDim,1))/dist;
                this.obs_matrix(obs_index, 2*DIM*(agent_id_i-1)+DIM+iDim) = 0;
                % Agent ID: j
                this.obs_matrix(obs_index, 2*DIM*(agent_id_j-1)+iDim)     = (pos_j(iDim,1)-pos_i(iDim,1))/dist;
                this.obs_matrix(obs_index, 2*DIM*(agent_id_j-1)+DIM+iDim) = 0;
            end
        end

        % Getters -------------------------------------------------------

        function output = getCurrentJointInformationVector(this)
            output = this.obs_info_vector_joint;
        end
        function output = getPreviousJointInformationVector(this)
            output = this.obs_info_vector_joint_prev;
        end
        function output = getCurrentJointInformationMatrix(this)
            output = this.obs_info_matrix_joint;
        end
        function output = getPreviousJointInformationMatrix(this)
            output = this.obs_info_matrix_joint_prev;
        end
    end
end