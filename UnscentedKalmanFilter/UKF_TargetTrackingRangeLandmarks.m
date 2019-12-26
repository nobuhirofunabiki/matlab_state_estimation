classdef UKF_TargetTrackingRangeLandmarks < UKF_LinearDynamics & RangeMeasurementLandmarks
    properties
    end

    methods
        function obj = UKF_TargetTrackingRangeLandmarks(args)
            obj@UKF_LinearDynamics(args.ukf);
            obj@RangeMeasurementLandmarks(args.rml);
        end

        function updateByMeasurement(this)
            num_sigma_points = size(this.sigma_points, 2);
            num_measures = size(this.getLandmarkList(), 2);
            z_pred = zeros(num_measures,1);
            Z = zeros(num_measures,num_sigma_points);
            obs_covmat = this.getMeasureCovarinaceMatrix();
            for iPoints = 1:num_sigma_points
                % TODO: only for 2D
                position = this.sigma_points(1:2, iPoints);
                this.setMeasurementVectorWithoutNoise(position);
                Z(:,iPoints) = this.getMeasurements();
                z_pred = z_pred + this.weights_mean(1,iPoints)*Z(:,iPoints);
            end
            this.z_pred = z_pred;
            this.Z_diff = Z - z_pred(:,ones(1,num_sigma_points));
            this.Pzz = this.Z_diff*diag(this.weights_cov)*(this.Z_diff)' + obs_covmat;
        end
    end
end