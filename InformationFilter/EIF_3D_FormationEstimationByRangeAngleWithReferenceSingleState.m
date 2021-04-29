classdef EIF_3D_FormationEstimationByRangeAngleWithReferenceSingleState < ...
    ExtendedInformationFilter

    properties (SetAccess = immutable)
        % Abstract properties of ExtendedInformationFilter
        process_noise_covmat
        discrete_system_matrix
        % Abstract properties of MultiagentUtilityBase
        num_agents
        num_dimensions
    end
    properties (SetAccess = protected)
        % Abstract properties of InformationFilterBase
        state_vector
        state_covmat
    end
    properties (SetAccess = private)
        range_sensor_       % Instance of RangeMeasurementMultiAgentWithReference class
        angle_sensor_       % Instance of AngleMeasurementMultiAgentWithReference3D class
        % Ratio actual noise vs noise model in the estimator
        ratio_range_noise
        ratio_angle_noise
        agent_id
    end

    methods (Access = public)
        function obj = EIF_3D_FormationEstimationByRangeAngleWithReferenceSingleState(args)
            obj@ExtendedInformationFilter(args);
            num_dims = args.num_dimensions;
            obj.range_sensor_           = RangeMeasurementMultiAgentWithReference(args.range_sensor);
            obj.angle_sensor_           = AngleMeasurementMultiAgentWithReference3D(args.angle_sensor);
            obj.ratio_range_noise       = args.ratio_range_noise;
            obj.ratio_angle_noise       = args.ratio_angle_noise;
            obj.num_agents              = args.num_agents;
            obj.num_dimensions          = num_dims;
            obj.state_vector            = args.state_vector;
            obj.agent_id                = args.agent_id;
            obj.process_noise_covmat    = args.process_noise_covmat;
            obj.discrete_system_matrix  = args.discrete_system_matrix;
            obj.state_covmat            = zeros(2*num_dims, 2*num_dims);
            for iDims = 1:num_dims
                obj.state_covmat(iDims, iDims) = args.sigma_position^2;
                obj.state_covmat(num_dims+iDims, num_dims+iDims) = args.sigma_velocity^2;
            end
        end
    end

    methods (Access = public)
        function executeInformationFilter(this, measures, discrete_system_matrix, adjacent_matrix, positions, position_ref)
            this.predictStateVectorAndCovariance(discrete_system_matrix);
            this.convertMomentsToInformationForm();

            num_dims    = this.num_dimensions;
            num_agents  = this.num_agents;

            for iAgents = 1:this.num_agents+1
                for jAgents = 1:this.num_agents+1
                    if (iAgents~=this.agent_id && jAgents~=this.agent_id)
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
            index_start = 2*num_dims*(this.agent_id-1)+1;
            index_end = index_start+2*num_dims-1;
            obs_matrix_range = obs_matrix_range_all(:,index_start:index_end);
            obs_covmat_range = this.range_sensor_.getMeasureCovarinaceMatrix()*this.ratio_range_noise;
            measures_predicted_range = this.range_sensor_.getMeasurements();
            this.addObservationInformation(...
                obs_matrix_range, obs_covmat_range, measures.ranges, measures_predicted_range);

            % Angle measurements
            this.angle_sensor_.computeMeasurementVector(positions, position_ref, false);
            this.angle_sensor_.setObservationMatrix(positions, position_ref);
            this.angle_sensor_.setMeasurementCovarianceMatrix(adjacent_matrix.angle);
            obs_matrix_angle_all = this.angle_sensor_.getObservationMatrix();
            index_start = 2*num_dims*(this.agent_id-1)+1;
            index_end = index_start+2*num_dims-1;
            obs_matrix_angle = obs_matrix_angle_all(:,index_start:index_end);
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
            this.addObservationInformation(...
                obs_matrix_angle, obs_covmat_angle, measures.angles, measures_predicted_angle);

            this.convertInformationToMomentsForm();
        end

        % Getters -------------------------------------------
        function output = getPositionVector(this)
            num_dims = this.num_dimensions;
            output = zeros(num_dims*this.num_agents,1);
            for iAgents = 1:this.num_agents
                output(num_dims*(iAgents-1)+1:num_dims*iAgents,1) ...
                    = this.state_vector(2*num_dims*(iAgents-1)+1:2*num_dims*(iAgents-1)+num_dims,1);
            end
        end
    end

    methods (Access = protected)
        % Override
        function predictStateVectorAndCovariance(this, discrete_system_matrix);
            Ad = discrete_system_matrix;
            Q = this.process_noise_covmat;
            this.state_vector = Ad*this.state_vector;
            this.state_covmat = Ad*this.state_covmat*Ad.' + Q;
        end
    end
    
end