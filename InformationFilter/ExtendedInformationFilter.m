classdef ExtendedInformationFilter < InformationFilter
    properties (SetAccess = protected)
        
    end

    methods
        function obj = ExtendedInformationFilter(args)
            obj@InformationFilter(args);
        end

        function executeInformationFilter(this, args)
            this.predictStateVectorAndCovariance();
            this.convertMomentsToInformationForm();
            this.addObservationInformation(...
                args.obs_matrix, args.obs_covmat, args.measures, args.measures_predicted);
            this.convertInformationToMomentsForm();
        end

        function addObservationInformation(this, ...
            obs_matrix, obs_covmat, measures, measures_predicted)
            H = obs_matrix;
            R = obs_covmat;
            y = measures;
            y_hat = measures_predicted;
            x_hat = this.state_vector;
            this.info_vector = this.info_vector + H.'/R*(y - y_hat + H*x_hat);
            this.info_matrix = this.info_matrix + H.'/R*H;
        end
    end

end