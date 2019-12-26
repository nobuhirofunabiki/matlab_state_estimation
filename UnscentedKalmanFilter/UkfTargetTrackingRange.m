classdef UkfTargetTrackingRange < UnscentedKalmanFilterBase
    properties (SetAccess = protected)
        range_list      % List of measurement source [row: position vector, column: source serial number]
        obs_cov
        discrete_system_matrix
        sys_cov
    end
    methods
        function obj = UkfTargetTrackingRange(args)
            obj@UnscentedKalmanFilterBase(args);
            obj.range_list = args.range_list;
            obj.obs_cov = args.obs_cov;
            obj.discrete_system_matrix = args.discrete_system_matrix;
            obj.sys_cov = args.sys_cov;
        end

        function predictStates(this)
            num_sigma_points = size(this.sigma_points, 2);
            state_est_prior = zeros(size(this.state_est));
            sigma_point_states = zeros(size(this.sigma_points));
            for iPoints = 1:num_sigma_points
                sigma_point_states(:,iPoints) = this.discrete_system_matrix*this.sigma_points(:,iPoints);
                state_est_prior = state_est_prior ...
                    + this.weights_mean(1,iPoints)*sigma_point_states(:,iPoints);
            end
            this.state_est = state_est_prior;
            this.X_diff = sigma_point_states - state_est_prior(:,ones(1,num_sigma_points));
            this.state_cov = this.X_diff*diag(this.weights_cov)*(this.X_diff)' + this.sys_cov;
        end

        function updateByMeasurement(this)
            num_sigma_points = size(this.sigma_points, 2);
            num_measures = size(this.range_list,2);
            z_pred = zeros(num_measures,1);
            Z = zeros(num_measures,num_sigma_points);
            for iPoints = 1:num_sigma_points
                % TODO: only for 2D
                position = this.sigma_points(1:2, iPoints);
                Z(:,iPoints) = this.calculateRangeMeasurements(position);
                z_pred = z_pred + this.weights_mean(1,iPoints)*Z(:,iPoints);
            end
            this.z_pred = z_pred;
            this.Z_diff = Z - z_pred(:,ones(1,num_sigma_points));
            this.Pzz = this.Z_diff*diag(this.weights_cov)*(this.Z_diff)' + this.obs_cov;
        end

        function output = calculateRangeMeasurements(this, position)
            num_measures = size(this.range_list,2);
            output = zeros(num_measures,1);
            for iMeasures = 1:num_measures
                output(iMeasures) = norm(position - this.range_list(:,iMeasures));
            end
        end

    end
end