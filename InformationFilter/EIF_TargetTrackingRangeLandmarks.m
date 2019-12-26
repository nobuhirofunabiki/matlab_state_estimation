classdef EIF_TargetTrackingRangeLandmarks < ExtendedInformationFilter & RangeMeasurementLandmarks
    properties (SetAccess = protected)
        
    end

    methods
        function obj = EIF_TargetTrackingRangeLandmarks(args)
            obj@ExtendedInformationFilter(args.eif);
            obj@RangeMeasurementLandmarks(args.rml);
        end

        function executeInformationFilter(this, measures)
            this.predictStateVectorAndCovariance();
            this.setMeasurementVectorWithoutNoise(this.getPosition());
            this.setObservationMatrix(this.getPosition());
            obs_matrix = this.getObservationMatrix();
            obs_covmat = this.getMeasureCovarinaceMatrix();
            measures_predicted = this.getMeasurements();
            this.convertMomentsToInformationForm();
            this.addObservationInformation(...
                obs_matrix, obs_covmat, measures, measures_predicted);
            this.convertInformationToMomentsForm();
        end

        % Getters --------------------------------------
        function output = getPosition(this)
            output = this.state_vector(1:this.num_dims, :);
        end
    end
end