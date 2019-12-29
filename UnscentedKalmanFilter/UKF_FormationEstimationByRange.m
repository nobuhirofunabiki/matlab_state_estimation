classdef UKF_FormationEstimationByRange < ...
    UKF_LinearDynamics & ...
    MultiagentUtilityBase

    properties (SetAccess = protected)
        state_vector
        state_covmat
        discrete_system_matrix
        process_noise_covmat
    end
    properties (SetAccess = private)
        range_sensor_       % Instance of RangeMeasurementInterAgents class
        position_sensor_    % Instance of PositionMeasurementMultiAgents class
    end
    properties (SetAccess = immutable)
        num_dimensions
        num_agents
    end

    methods (Access = public)
        function obj = UKF_FormationEstimationByRange(args)
            obj@UKF_LinearDynamics(args);
            obj.range_sensor_           = RangeMeasurementInterAgents(args.rmia);
            obj.position_sensor_        = PositionMeasurementMultiAgents(args.pmb);
            obj.state_vector            = args.init_state_vector;
            obj.state_covmat            = obj.createStateCovarianceMatrix(args.sigma_position, args.sigma_velocity);
            obj.process_noise_covmat    = args.process_noise_covmat;
            obj.discrete_system_matrix  = obj.createDiscreteSystemMatrix(args.discrete_system_matrix);
            obj.num_dimensions          = args.num_dimensions;
            obj.num_agents              = args.num_agents;
        end
    end

    methods (Access = public)
        function executeUnscentedKalmanFilter(this, measurements, adjacent_matrix)
            this.generateSigmaPoints();
            this.predictStates();
            this.generateSigmaPoints();
            this.processMeasurements(measurements, adjacent_matrix);
            this.updateByMeasurements();
        end
    end

    methods (Access = protected)
        function processMeasurements(this, measurements, adjacent_matrix)
            this.measurements   = vertcat(measurements.ranges, measurements.positions);
            num_sigma_points    = size(this.sigma_points, 2);
            num_measures        = size(this.measurements, 1);
            num_ranges          = size(measurements.ranges, 1);
            num_positions       = size(measurements.positions, 1);
            z_pred              = zeros(num_measures,1);
            Z                   = zeros(num_measures,num_sigma_points);
            obs_covmat          = zeros(num_measures, num_measures);
            this.range_sensor_.updateMeasurementCovarianceMatrix(adjacent_matrix);
            obs_covmat(1:num_ranges, 1:num_ranges) ...
                = this.range_sensor_.getMeasureCovarinaceMatrix();
            obs_covmat(num_ranges+1:num_measures, num_ranges+1:num_measures) ...
                = this.position_sensor_.getMeasureCovarinaceMatrix();
            for iPoints = 1:num_sigma_points
                positions = this.getPositionVector(iPoints);
                this.range_sensor_.setMeasurementVectorWithoutNoise(positions);
                this.position_sensor_.computeMeasurementVector(positions, false);
                Z(1:num_ranges,iPoints)         = this.range_sensor_.getMeasurements();
                Z(num_ranges+1:num_measures)    = this.position_sensor_.getMeasurements();
                z_pred = z_pred + this.weights_mean(1,iPoints)*Z(:,iPoints);
            end
            this.z_pred = z_pred;
            this.Z_diff = Z - z_pred(:,ones(1,num_sigma_points));
            this.Pzz    = this.Z_diff*diag(this.weights_cov)*(this.Z_diff)' + obs_covmat;
        end
    end

    methods (Access = private)
        function output = getPositionVector(this, iSigmaPoints)
            num_dims = this.num_dimensions;
            output = zeros(num_dims*this.num_agents,1);
            for iAgents = 1:this.num_agents
                output(num_dims*(iAgents-1)+1:num_dims*iAgents,1) ...
                    = this.sigma_points(2*num_dims*(iAgents-1)+1:2*num_dims*(iAgents-1)+num_dims,iSigmaPoints);
            end
        end
    end

end