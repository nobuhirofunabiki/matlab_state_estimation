classdef ExtendedInformationFilter < InformationFilterBase
    properties (Abstract = true, SetAccess = immutable)
        process_noise_covmat
        discrete_system_matrix
    end

    methods (Access = protected)
        function obj = ExtendedInformationFilter(args)
            obj@InformationFilterBase(args);
        end
    end

    methods (Access = protected)
        function predictStateVectorAndCovariance(this)
            Ad = this.discrete_system_matrix;
            Q = this.process_noise_covmat;
            this.state_vector = Ad*this.state_vector;
            this.state_covmat = Ad*this.state_covmat*Ad.' + Q;
        end

        function addObservationInformation(this, ...
            obs_matrix, obs_covmat, measures, measures_predicted)
            H = obs_matrix;
            R = obs_covmat;
            y = measures;
            y_hat = measures_predicted;
            x_hat = this.state_vector;
            this.obs_info_vector    = this.obs_info_vector + H.'/R*(y - y_hat + H*x_hat);
            this.obs_info_matrix    = this.obs_info_vector + H.'/R*H;
            this.info_vector        = this.info_vector + H.'/R*(y - y_hat + H*x_hat);
            this.info_matrix        = this.info_matrix + H.'/R*H;
        end
    end

end