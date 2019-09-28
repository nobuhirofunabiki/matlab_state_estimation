classdef RangeKalmanFilter < ExtendedKalmanFilter
    properties (SetAccess = private)
        
    end
    methods
        function obj = RangeKalmanFilter(args)
            obj@ExtendedKalmanFilter(args);
        end

        % Set functions
        function setFilterParameters(this, agent_, dynamics_, measurement_, ...
                                     num_edges, init_error_sigma)
            NUM_AGENTS = size(agent_, 2);
            Ad = dynamics_.getDiscreteSystemMatrix;
            num_var = this.num_variable;
            num_pos = num_var/2;
            num_measure = this.num_measure;

            A = zeros(6*NUM_AGENTS, 6*NUM_AGENTS);
            for i = 1:NUM_AGENTS
                A(1+6*(i-1):6*i,1+6*(i-1):6*i) = Ad;
            end

            P = zeros(num_var, num_var);
            for i = 1:num_pos
                P(i, i) = (init_error_sigma.position)^2;
                P(num_pos+i, num_pos+i) = (init_error_sigma.velocity)^2;
            end

            SIGMA_RANGE = measurement_.getNoiseSigma();
            R = diag(SIGMA_RANGE^2*ones(num_edges, 1));

            this.setDiscreteSystemMatrix(A);
            this.setStateCovarianceMatrix(P);
            this.setMeasurementCovarianceMatrix(R);
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
        function setObservationMatrixReference(this, pos_chief, pos_i, agent_id_i, obs_index)
            dist = norm(pos_i - pos_chief);
            DIM  = numel(pos_i);

            for iDim = 1:DIM
                this.obs_matrix(obs_index, 2*DIM*(agent_id_i-1)+iDim)     = (pos_i(iDim,1)-pos_chief(iDim,1))/dist;
                this.obs_matrix(obs_index, 2*DIM*(agent_id_i-1)+DIM+iDim) = 0;
            end
        end

    end
end