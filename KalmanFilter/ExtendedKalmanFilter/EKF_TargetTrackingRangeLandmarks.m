classdef EKF_TargetTrackingRangeLandmarks < ExtendedKalmanFilterBase
    properties (SetAccess = protected)
        range_sensor_;      % Class instance of RangeMeasurementLandmarks
        state_vector
        state_covmat
        process_noise_covmat
        discrete_system_matrix
    end

    properties (SetAccess = immutable)
        num_dimensions
    end

    methods (Access = public)
        function obj = EKF_TargetTrackingRangeLandmarks(args)
            obj@ExtendedKalmanFilterBase(args);
            obj.range_sensor_           = RangeMeasurementLandmarks(args.rml);
            obj.state_vector            = args.state_vector;
            obj.state_covmat            = args.state_covmat;
            obj.process_noise_covmat    = args.process_noise_covmat;
            obj.discrete_system_matrix  = args.discrete_system_matrix;
            obj.num_dimensions          = args.num_dimensions;
        end
    end

    methods (Access = protected)
        function processMeasurements(this, measurements)
            this.range_sensor_.computeMeasurementVector(this.getPosition(), false);
            this.range_sensor_.setObservationMatrix(this.getPosition());
            H = this.range_sensor_.getObservationMatrix();
            R = this.range_sensor_.getMeasureCovarinaceMatrix();
            z_hat = this.range_sensor_.getMeasurements();
            z = measurements;
            x_hat = this.state_vector;
            P = this.state_covmat;
            K = P*H.'/(H*P*H.'+R);
            this.state_vector = x_hat + K*(z-z_hat);
            this.state_covmat = (eye(size(P))-K*H)*P;
        end
    end

    methods (Access = private)
        % Getters --------------------------------------
        function output = getPosition(this)
            output = this.state_vector(1:this.num_dimensions, :);
        end
    end
end