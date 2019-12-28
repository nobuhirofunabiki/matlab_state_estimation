classdef DIF_FormationEstimationByRange < DIF_LinearDynamics
    properties (SetAccess = private)
        range_sensor_       % Instance of RangeMeasurementInterAgents class
        position_sensor_    % Instance of PositionMeasurementMultiAgents class
    end

    methods (Access = public)
        function obj = DIF_FormationEstimationByRange(args)
            obj@DIF_LinearDynamics(args);
            obj.range_sensor_       = RangeMeasurementInterAgents(args.rmia);
            obj.position_sensor_    = PositionMeasurementMultiAgents(args.pmb);
        end
    end

    methods (Access = protected)
        function processMeasurements(this, measures, adjacent_matrix)
            % measures should have 'ranges' and 'positions' field
            positions = this.getPositionVector();

            % Range measurements
            this.range_sensor_.calculateMeasurementVector(positions);
            this.range_sensor_.setObservationMatrix(positions);
            this.range_sensor_.updateMeasurementCovarianceMatrix(adjacent_matrix);
            obs_matrix_range = this.range_sensor_.getObservationMatrix();
            obs_covmat_range = this.range_sensor_.getMeasureCovarinaceMatrix();
            measures_predicted_range = this.range_sensor_.getMeasurements();
            this.addObservationInformationRange(...
                obs_matrix_range, obs_covmat_range, measures.ranges, measures_predicted_range);
            
            % Position measurements
            this.position_sensor_.computeMeasurementVector(positions, false);
            obs_matrix_pos = this.position_sensor_.getObservationMatrix;
            obs_covmat_pos = this.position_sensor_.getMeasureCovarinaceMatrix;
            % measures_predicted_pos = this.position_sensor_.getMeasurements();
            this.addObservationInformation(...
                obs_matrix_pos, obs_covmat_pos, measures.positions);   
        end
    end

    methods(Access = private)
        function addObservationInformationRange(this,...
            obs_matrix, obs_covmat, measures, measures_predicted)
            H = obs_matrix;
            R = obs_covmat;
            y = measures;
            y_hat = measures_predicted;
            x_hat = this.state_vector;
            this.obs_info_vector = this.obs_info_vector + H.'/R*(y - y_hat + H*x_hat);
            this.obs_info_matrix = this.obs_info_matrix + H.'/R*H;
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

end