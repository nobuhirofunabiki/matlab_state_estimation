classdef ExtendedKalmanFilterBase < handle
    properties (Abstract = true, SetAccess = protected)
        state_vector
        state_covmat
        process_noise_covmat
        discrete_system_matrix
    end

    methods(Access = protected)
        function obj = ExtendedKalmanFilterBase(args)
        end

        function predictStateVectorAndCovariance(this)
            Ad = this.discrete_system_matrix;
            Q = this.process_noise_covmat;
            this.state_vector = Ad*this.state_vector;
            this.state_covmat = Ad*this.state_covmat*Ad.' + Q;
        end
    end

    methods (Abstract = true, Access = protected)
        processMeasurements(this);
    end

    methods (Access = public)
        function executeFiltering(this, measurements)
            this.predictStateVectorAndCovariance();
            this.processMeasurements(measurements);
        end

        % Setters-------------------------------------------------------
        function setEstimatedVariable(this, index_begin, index_end, arg_variable)
            this.state_vector(index_begin:index_end, 1) = arg_variable;
        end

        function setStateVectorAll(this, arg)
            this.state_vector = arg;
        end

        function setStateCovarianceMatrixAll(this, state_covmat)
            this.state_covmat = state_covmat;
        end

        % Getters -------------------------------------------------------
        function output = getStateVector(this)
            output = this.state_vector;
        end

        function output = getStateCovarianceMatrix(this)
            output = this.state_covmat;
        end
    end
end