classdef UKF_TargetTrackingRangeLandmarks < UKF_LinearDynamics
    properties (SetAccess = protected)
        state_vector
        state_covmat
        discrete_system_matrix
        process_noise_covmat 
    end
    properties (SetAccess = private)
        range_sensor_       % Class instance of RangeMeasurementLandmarks
    end

    methods (Access = public)
        function obj = UKF_TargetTrackingRangeLandmarks(args)
            obj@UKF_LinearDynamics(args);
            obj.range_sensor_           = RangeMeasurementLandmarks(args.rml);
            obj.state_vector            = args.state_vector;
            obj.state_covmat            = args.state_covmat;
            obj.discrete_system_matrix  = args.discrete_system_matrix;
            obj.process_noise_covmat    = args.process_noise_covmat;
        end
    end

    methods (Access = protected)
        function processMeasurements(this, measurements)
            this.measurements   = measurements;
            num_sigma_points    = size(this.sigma_points, 2);
            num_measures        = size(this.range_sensor_.getLandmarkList(), 2);
            z_pred              = zeros(num_measures,1);
            Z                   = zeros(num_measures,num_sigma_points);
            obs_covmat          = this.range_sensor_.getMeasureCovarinaceMatrix();
            for iPoints = 1:num_sigma_points
                % TODO: only for 2D
                position = this.sigma_points(1:2, iPoints);
                this.range_sensor_.computeMeasurementVector(position, false);
                Z(:,iPoints) = this.range_sensor_.getMeasurements();
                z_pred = z_pred + this.weights_mean(1,iPoints)*Z(:,iPoints);
            end
            this.z_pred = z_pred;
            this.Z_diff = Z - z_pred(:,ones(1,num_sigma_points));
            this.Pzz    = this.Z_diff*diag(this.weights_cov)*(this.Z_diff)' + obs_covmat;
        end
    end
end