classdef EIF_TargetTrackingRangeLandmarks < ExtendedInformationFilter
    properties (SetAccess = protected)
        range_sensor_;      % Class instance of RangeMeasurementLandmarks
        state_vector        % x(k): State vector
        state_covmat        % P(k): Covariance matrix of state
    end
    properties (SetAccess = immutable)
        process_noise_covmat
        discrete_system_matrix
        num_dimensions
    end

    methods (Access = public)
        function obj = EIF_TargetTrackingRangeLandmarks(args)
            obj@ExtendedInformationFilter(args);
            % obj.checkConstructorArguments(args);
            obj.range_sensor_           = RangeMeasurementLandmarks(args.rml);
            obj.state_vector            = args.state_vector;
            obj.state_covmat            = args.state_covmat;
            obj.process_noise_covmat    = args.process_noise_covmat;
            obj.discrete_system_matrix  = args.discrete_system_matrix;
            obj.num_dimensions          = args.num_dimensions;
        end

        function executeInformationFilter(this, measures)
            this.predictStateVectorAndCovariance();
            this.range_sensor_.computeMeasurementVector(this.getPosition(), false);
            this.range_sensor_.setObservationMatrix(this.getPosition());
            obs_matrix = this.range_sensor_.getObservationMatrix();
            obs_covmat = this.range_sensor_.getMeasureCovarinaceMatrix();
            measures_predicted = this.range_sensor_.getMeasurements();
            this.convertMomentsToInformationForm();
            this.addObservationInformation(...
                obs_matrix, obs_covmat, measures, measures_predicted);
            this.convertInformationToMomentsForm();
        end
    end

    methods (Access = private)
        function checkConstructorArguments(this, args)
            disp("Check constructor arguments for EIF_TargetTrackingRangeLandmarks");
            num_vars = this.num_variables;
            assert(isequal(size(args.state_vector), [num_vars, 1]), ...
                "state_vector is NOT a correct size vector");
            assert(isequal(size(args.state_covmat), [num_vars, num_vars]), ...
                "state_covmat is NOT a correct size matrix");
            assert(isequal(size(args.process_noise_covmat), [num_vars, num_vars]), ...
                "process_noise_covmat is NOT a correct size matrix");
            assert(isequal(size(args.discrete_system_matrix), [num_vars, num_vars]), ...
                "discrete_system_matrix is NOT a correct size matrix");
        end

        % Getters --------------------------------------
        function output = getPosition(this)
            output = this.state_vector(1:this.num_dimensions, :);
        end
    end
end