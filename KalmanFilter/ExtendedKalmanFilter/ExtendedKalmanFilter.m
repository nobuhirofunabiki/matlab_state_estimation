classdef ExtendedKalmanFilter < handle
    properties (SetAccess = protected)
        num_variable;
        num_measure;
        estimated_variable;
        obs_matrix;
        obs_measured;
        obs_estimated;
        covmat_state;
        covmat_measure;
        Ad;
    end
    methods
        function obj = ExtendedKalmanFilter(args)
            num_variable = args.num_variable;
            num_measure  = args.num_measure;
            obj.num_variable        = num_variable;
            obj.num_measure         = num_measure;
            obj.estimated_variable  = zeros(num_variable, 1);
            obj.obs_matrix          = zeros(num_measure, num_variable);
            obj.obs_measured        = zeros(num_measure, 1);
            obj.obs_estimated       = zeros(num_measure, 1);
            obj.covmat_state        = zeros(num_variable, num_variable);
            obj.covmat_measure      = zeros(num_measure, num_measure);
            obj.Ad                  = zeros(num_variable, num_variable);
        end
        function x_hat = executeEKF(this)
            x_pri   = this.estimated_variable;
            P_pri   = this.covmat_state;
            C       = this.obs_matrix;
            R       = this.covmat_measure;
            D       = this.obs_measured;
            D_est   = this.obs_estimated;
            G = P_pri*C.'/(C*P_pri*C.'+R);
            x_hat = x_pri+G*(D-D_est);
            this.covmat_state = (eye(size(P_pri))-G*C)*P_pri;
        end
        function propagateCovarianceMatrix(this)
            Ad = this.Ad;
            % TODO: Add control and disturbance noises
            this.covmat_state = Ad*this.covmat_state*Ad.';
        end
        function setEstimatedVariable(this, index_begin, index_end, arg_variable)
            if length(arg_variable) ~= (index_end-index_begin)+1
                error('Index designation is invaild.')
            end
            this.estimated_variable(index_begin:index_end, 1) = arg_variable;
        end
        function setStateCovarianceMatrix(this, arg_covmat_state)
            this.covmat_state = arg_covmat_state;
        end
        function setMeasurementData(this, arg_obs_measured, arg_obs_estimated)
            this.obs_measured  = arg_obs_measured;
            this.obs_estimated = arg_obs_estimated;
        end
        function setDiscreteSystemMatrix(this, arg_system_matrix)
            this.Ad = arg_system_matrix;
        end
        function setMeasurementCovarianceMatrix(this, arg_covmat_measure)
            this.covmat_measure = arg_covmat_measure;
        end
    end
end