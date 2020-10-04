classdef DecentralizedSquareRootUnscentedInformationFilterBase < DecentralizedUnscentedInformationFilterBase
    properties (Abstract = true, SetAccess = protected)
        S
        Sz
        sqrtQ
        sqrtR
    end

    properties (SetAccess = protected)
        functions_ = CommonEstimatorFunctions()
    end

    methods (Access = protected)
        function obj = DecentralizedSquareRootUnscentedInformationFilterBase(args)
            obj@DecentralizedUnscentedInformationFilterBase(args);
        end

        function output = generateSigmaPoints(this)
            A = this.gamma * (this.S).';
            x = this.state_vector;
            Y = x(:,ones(1,numel(x)));
            output = [x Y+A Y-A];
        end

        function predictStateVectorAndCovariance(this)
            num_sigma_points = size(this.sigma_points, 2);
            state_vector_prior = zeros(size(this.state_vector));
            sigma_point_states = zeros(size(this.sigma_points));
            for iPoints = 1:num_sigma_points
                sigma_point_states(:,iPoints) = this.discrete_system_matrix*this.sigma_points(:,iPoints);
                state_vector_prior = state_vector_prior ...
                    + this.weights_mean(1,iPoints)*sigma_point_states(:,iPoints);
            end
            this.state_vector = state_vector_prior;
            this.X_diff = sigma_point_states - state_vector_prior(:,ones(1,num_sigma_points));
            this.S = this.functions_.qr_decomposition(this.X_diff, this.sqrtQ, this.weights_cov);
        end
    end
end