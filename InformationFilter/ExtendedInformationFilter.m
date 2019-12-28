classdef ExtendedInformationFilter < InformationFilterBase
    properties (SetAccess = protected)
        
    end

    methods (Access = protected)
        function obj = ExtendedInformationFilter(args)
            obj@InformationFilterBase(args);
        end
    end

    methods (Access = protected)
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