classdef ExtendedInformationFilter < InformationFilter
    properties (SetAccess = protected)
        
    end

    methods
        function obj = ExtendedInformationFilter(args)
            obj@InformationFilter(args);
        end

        function addObservationInformation(this, args)
            H = args.obs_matrix;
            R = args.obs_covmat;
            y = args.measures;
            y_hat = args.measures_predicted;
            x_hat = this.state_vector;
            this.info_vector = this.info_vector + H.'/R*(y - y_hat + H*x_hat);
            this.info_matrix = this.info_matrix + H.'/R*H;
        end
    end

end