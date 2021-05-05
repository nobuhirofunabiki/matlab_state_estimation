classdef PartialStateEstimation < MultiagentUtilityBase

    properties (SetAccess = immutable)
        num_agents
        num_dimensions
    end
    properties (SetAccess = protected)
        state_vector
        state_covmat
        current_estimated_agents
        previous_estimated_agents
        num_estimated
        agent_id
        process_noise_covmat
    end
    properties (SetAccess = private)
        range_sensor_       % Instance of RangeMeasurementMultiAgentWithReference class
        angle_sensor_       % Instance of AngleMeasurementMultiAgentWithReference3D class
        % Ratio actual noise vs noise model in the estimator
        ratio_range_noise
        ratio_angle_noise
    end

    methods (Access = public)
        function obj = PartialStateEstimation(args)
            num_agents                      = args.num_agents;
            obj.num_agents                  = num_agents;
            obj.num_dimensions              = args.num_dimensions;
            obj.current_estimated_agents    = zeros(1,num_agents);
            obj.previous_estimated_agents   = zeros(1,num_agents);
            obj.agent_id                    = args.agent_id;
            obj.process_noise_covmat        = args.process_noise_covmat;
            obj.state_vector                = args.state_vector;
            obj.state_covmat                = obj.createStateCovarianceMatrix(args.sigma_position, args.sigma_velocity);
            obj.range_sensor_               = RangeMeasurementMultiAgentWithReference(args.range_sensor);
            obj.angle_sensor_               = AngleMeasurementMultiAgentWithReference3D(args.angle_sensor);
            obj.ratio_range_noise           = args.ratio_range_noise;
            obj.ratio_angle_noise           = args.ratio_angle_noise;
        end
    end

    methods (Access = public)
        function executeInformationFilter(this, measurements, discrete_system_matrix, adjacent_matrix, position_ref, shared_state_vector, shared_state_covmat)
            this.checkNetworkUpdate(adjacent_matrix.comm);
            state_vector = this.makeStateVector(shared_state_vector);
            state_covmat = this.makeStateCovarianceMatrix(shared_state_covmat);
            [state_vector, state_covmat] = this.predictStateVectorAndCovariance(discrete_system_matrix, state_vector, state_covmat);
            [info_vector, info_matrix] = this.convertMomentsToInformationForm(state_vector, state_covmat);

            num_dims = this.num_dimensions;
            positions = zeros(num_dims*this.num_agents,1);
            for iAgents = 1:this.num_agents
                positions(num_dims*(iAgents-1)+1:num_dims*iAgents,1) ...
                    = shared_state_vector(2*num_dims*(iAgents-1)+1:2*num_dims*(iAgents-1)+num_dims,1);
            end

            [info_vector, info_matrix] = this.processMeasurementData(state_vector, info_vector, info_matrix, measurements, adjacent_matrix, positions, position_ref);
            [state_vector, state_covmat] = this.convertInformationToMomentsForm(info_vector, info_matrix);
            this.updateStateVector(state_vector);
            this.updateStateCovarianceMatrix(state_covmat);
        end
    end

    methods (Access = protected)
        function checkNetworkUpdate(this, adjacent_matrix)
            % Need to remove the reference spacecraft from the adjacent matrix
            query_estimated_agents = adjacent_matrix(this.agent_id, 1:this.num_agents);
            query_estimated_agents(1,this.agent_id) = 1;
            if not(isequal(query_estimated_agents, this.current_estimated_agents))
                this.previous_estimated_agents = this.current_estimated_agents;
                this.current_estimated_agents = query_estimated_agents;
                this.num_estimated = sum(this.current_estimated_agents);
            end
        end

        function [state_vector, state_covmat] = predictStateVectorAndCovariance(this, discrete_system_matrix_for_one_agent, state_vector, state_covmat)
            num_dims = this.num_dimensions;
            num_estimated = this.num_estimated;
            Ad = zeros(2*num_dims*num_estimated, 2*num_dims*num_estimated);
            Q = zeros(2*num_dims*num_estimated, 2*num_dims*num_estimated);
            for iAgents = 1:num_estimated
                Ad(1+2*num_dims*(iAgents-1):2*num_dims*iAgents,...
                    1+2*num_dims*(iAgents-1):2*num_dims*iAgents)...
                = discrete_system_matrix_for_one_agent;
                Q(1+2*num_dims*(iAgents-1):2*num_dims*iAgents,...
                    1+2*num_dims*(iAgents-1):2*num_dims*iAgents)...
                = this.process_noise_covmat;
            end
            state_vector = Ad*state_vector;
            state_covmat = Ad*state_covmat*Ad.' + Q;
        end

        function [info_vector, info_matrix] = addObservationInformation(this, ...
            state_vector, info_vector, info_matrix, obs_matrix, obs_covmat, measures, measures_predicted)
            H = obs_matrix;
            R = obs_covmat;
            y = measures;
            y_hat = measures_predicted;
            x_hat = state_vector;
            info_vector = info_vector + H.'/R*(y - y_hat + H*x_hat);
            info_matrix = info_matrix + H.'/R*H;
        end

        function [info_vector, info_matrix] = processMeasurementData(this, ...
            state_vector, info_vector, info_matrix, measures, adjacent_matrix, positions, position_ref)

            num_dims = this.num_dimensions;
            num_agents = this.num_agents;

            % [TODO] Need to consider how to dear with the reference spacecraft
            for iAgents = 1:num_agents
                for jAgents = 1:num_agents
                    b_iAgents_estimated = this.current_estimated_agents(1,iAgents);
                    b_jAgents_estimated = this.current_estimated_agents(1,jAgents);
                    if (b_iAgents_estimated == 0 && b_jAgents_estimated == 0)
                        adjacent_matrix.range(iAgents, jAgents) = 0;
                        adjacent_matrix.angle(iAgents, jAgents) = 0;
                    end
                end
            end

            % Range measurements
            this.range_sensor_.computeMeasurementVector(positions, position_ref, false);
            this.range_sensor_.setObservationMatrix(positions, position_ref);
            this.range_sensor_.updateMeasurementCovarianceMatrix(adjacent_matrix.range);
            obs_matrix_range_all = this.range_sensor_.getObservationMatrix();
            size_row = size(obs_matrix_range_all,1);
            size_column = 2*num_dims*this.num_estimated;
            obs_matrix_range = zeros(size_row, size_column);
            iEstimated = 0;
            for iAgents = 1:num_agents
                if this.current_estimated_agents(1,iAgents) == 1
                    iEstimated = iEstimated + 1;
                    index_start = 2*num_dims*(iEstimated-1)+1;
                    index_end = index_start+2*num_dims-1;
                    obs_matrix_range(:,index_start:index_end) ...
                        = obs_matrix_range_all(:,2*num_dims*(iAgents-1)+1:2*num_dims*iAgents);
                end
            end
            obs_covmat_range = this.range_sensor_.getMeasureCovarinaceMatrix()*this.ratio_range_noise;
            measures_predicted_range = this.range_sensor_.getMeasurements();
            [info_vector, info_matrix] = this.addObservationInformation(...
                state_vector, info_vector, info_matrix, obs_matrix_range, obs_covmat_range, measures.ranges, measures_predicted_range);
            
            % Angle measurements
            this.angle_sensor_.computeMeasurementVector(positions, position_ref, false);
            this.angle_sensor_.setObservationMatrix(positions, position_ref);
            this.angle_sensor_.setMeasurementCovarianceMatrix(adjacent_matrix.angle);
            obs_matrix_angle_all = this.angle_sensor_.getObservationMatrix();
            size_row = size(obs_matrix_angle_all,1);
            size_column = 2*num_dims*this.num_estimated;
            obs_matrix_angle = zeros(size_row, size_column);
            iEstimated = 0;
            for iAgents = 1:num_agents
                if this.current_estimated_agents(1,iAgents) == 1
                    iEstimated = iEstimated + 1;
                    index_start = 2*num_dims*(iEstimated-1)+1;
                    index_end = index_start+2*num_dims-1;
                    obs_matrix_angle(:,index_start:index_end) ...
                        = obs_matrix_angle_all(:,2*num_dims*(iAgents-1)+1:2*num_dims*iAgents);
                end
            end
            obs_covmat_angle = this.angle_sensor_.getMeasureCovarinaceMatrix()*this.ratio_angle_noise;
            measures_predicted_angle = this.angle_sensor_.getMeasurements();
            % TODO: How to tuckle the following singular point problem?
            diff_measures = measures_predicted_angle - measures.angles;
            for iMeasures = 1:length(diff_measures)
                if (abs(diff_measures(iMeasures,1)) >= pi)
                    if (measures_predicted_angle(iMeasures,1) > measures.angles(iMeasures,1))
                        measures_predicted_angle(iMeasures,1) = measures_predicted_angle(iMeasures,1) - 2*pi;
                    else
                        measures_predicted_angle(iMeasures,1) = measures_predicted_angle(iMeasures,1) + 2*pi;
                    end
                end
            end
            [info_vector, info_matrix] = this.addObservationInformation(...
                state_vector, info_vector, info_matrix, obs_matrix_angle, obs_covmat_angle, measures.angles, measures_predicted_angle);
        end

        function state_vector = makeStateVector(this, shared_state_vector)
            num_est = this.num_estimated;
            num_dims = this.num_dimensions;
            state_vector = zeros(2*num_dims*num_est, 1);
            iEstimated = 0;
            for iAgents = 1:this.num_agents
                if this.current_estimated_agents(1,iAgents) == 1
                    iEstimated = iEstimated + 1;
                    state_vector_one_agent = zeros(2*num_dims, 1);
                    if this.previous_estimated_agents(1,iAgents) == 1
                        state_vector_one_agent = this.state_vector(2*num_dims*(iAgents-1)+1:2*num_dims*iAgents, 1);
                    else
                        state_vector_one_agent = shared_state_vector(2*num_dims*(iAgents-1)+1:2*num_dims*iAgents, 1);
                    end
                    state_vector(2*num_dims*(iEstimated-1)+1:2*num_dims*iEstimated, 1) = state_vector_one_agent;
                end
            end
            if iEstimated ~= this.num_estimated
                error('Invalid counting number for makeStateVector');
                exit
            end
        end

        function state_covmat = makeStateCovarianceMatrix(this, shared_state_covmat)
            num_est = this.num_estimated;
            num_dims = this.num_dimensions;
            state_covmat = zeros(2*num_dims*num_est, 2*num_dims*num_est);
            iEstimated = 0;
            for iAgents = 1:this.num_agents
                if this.current_estimated_agents(1,iAgents) == 1
                    iEstimated = iEstimated + 1;
                    jEstimated = 0;
                    for jAgents = 1:this.num_agents
                        if this.current_estimated_agents(1,jAgents) == 1
                            jEstimated = jEstimated + 1;
                            state_covmat_one_agent = zeros(2*num_dims, 2*num_dims);
                            if this.previous_estimated_agents(1,iAgents) == 1 && this.previous_estimated_agents(1,jAgents) == 1
                                iAgents_pre_index = sum(this.previous_estimated_agents(1,1:iAgents));
                                jAgents_pre_index = sum(this.previous_estimated_agents(1,1:jAgents));
                                state_covmat_one_agent = this.state_covmat(2*num_dims*(iAgents_pre_index-1)+1:2*num_dims*iAgents_pre_index, 2*num_dims*(jAgents_pre_index-1)+1:2*num_dims*jAgents_pre_index);
                            elseif this.previous_estimated_agents(1,iAgents) ~= 1 && iAgents == jAgents
                                state_covmat_one_agent = shared_state_covmat(2*num_dims*(iAgents-1)+1:2*num_dims*iAgents, 2*num_dims*(jAgents-1)+1:2*num_dims*jAgents);
                            else
                                state_covmat_one_agent = zeros(2*num_dims, 2*num_dims);
                            end
                            state_covmat(2*num_dims*(iEstimated-1)+1:2*num_dims*iEstimated, 2*num_dims*(jEstimated-1)+1:2*num_dims*jEstimated) ...
                                = state_covmat_one_agent;
                        end
                    end
                end
            end
        end

        function [info_vector, info_matrix] = convertMomentsToInformationForm(this, state_vector, state_covmat)
            info_matrix = inv(state_covmat);
            info_vector =  info_matrix * state_vector;
        end

        function [state_vector, state_covmat] = convertInformationToMomentsForm(this, info_vector, info_matrix)
            state_covmat = inv(info_matrix);
            state_vector = state_covmat * info_vector;
        end

        function updateStateVector(this, state_vector)
            num_dims = this.num_dimensions;
            iEstimated = 0;
            for iAgents = 1:this.num_agents
                if this.current_estimated_agents(1,iAgents) == 1
                    iEstimated = iEstimated + 1;
                    this.state_vector(2*num_dims*(iAgents-1)+1:2*num_dims*iAgents, 1)...
                        = state_vector(2*num_dims*(iEstimated-1)+1:2*num_dims*iEstimated, 1);
                end
            end
        end

        function updateStateCovarianceMatrix(this, state_covmat)
            num_dims = this.num_dimensions;
            for iAgents = 1:this.num_agents
                for jAgents = 1:this.num_agents
                    b_iAgents_estimated = (this.current_estimated_agents(1,iAgents) == 1);
                    b_jAgents_estimated = (this.current_estimated_agents(1,jAgents) == 1);
                    if (b_iAgents_estimated==true && b_jAgents_estimated==true)
                        iEstimated = sum(this.current_estimated_agents(1,1:iAgents));
                        jEstimated = sum(this.current_estimated_agents(1,1:jAgents));
                        this.state_covmat(2*num_dims*(iAgents-1)+1:2*num_dims*iAgents, 2*num_dims*(jAgents-1)+1:2*num_dims*jAgents) ...
                            = state_covmat(2*num_dims*(iEstimated-1)+1:2*num_dims*iEstimated, 2*num_dims*(jEstimated-1)+1:2*num_dims*jEstimated);
                    end
                end
            end
        end
    end

    methods (Access = public)

        %%%%%%%%%%%%%%%%%%%
        % Setters
        %%%%%%%%%%%%%%%%%%%

        function setStateVector(this, arg_state_vector)
            this.state_vector = arg_state_vector;
        end

        function setStateCovarianceMatrix(this, arg_state_covmat)
            this.state_covmat = arg_state_covmat;
        end

        function setPreviousEstimatedAgents(this, arg)
            this.previous_estimated_agents = arg;
        end

        function setCurrentEstimatedAgents(this, arg)
            this.current_estimated_agents = arg;
        end

        function setNumberEstimated(this, arg)
            this.num_estimated = arg;
        end

        %%%%%%%%%%%%%%%%%%%
        % Getters
        %%%%%%%%%%%%%%%%%%%

        function output = getPreviousEstimatedAgents(this)
            output = this.previous_estimated_agents;
        end

        function output = getCurrentEstimatedAgents(this)
            output = this.current_estimated_agents;
        end

        function output = getNumberEstimated(this)
            output = this.num_estimated;
        end

        function output = getAgentID(this)
            output = this.agent_id;
        end

        function output = getStateVector(this)
            output = this.state_vector;
        end

        function output = getStateVectorOfAgentID(this)
            iAgents = this.agent_id;
            num_dims = this.num_dimensions;
            output = this.state_vector(2*num_dims*(iAgents-1)+1:2*num_dims*iAgents,1);
        end

        function output = getStateCovarianceMatrixOfAgentID(this)
            iAgents = this.agent_id;
            num_dims = this.num_dimensions;
            output = this.state_covmat(...
                2*num_dims*(iAgents-1)+1:2*num_dims*iAgents, 2*num_dims*(iAgents-1)+1:2*num_dims*iAgents);
        end
    end
end