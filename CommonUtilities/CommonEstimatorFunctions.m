classdef CommonEstimatorFunctions < handle

    methods (Access = public)
        function obj = CommonEstimatorFunctions(args)
        end

        function output = qr_decomposition(this, x_diff, square_root_covariance, weights_cov)
            num_sigma_points = size(x_diff,2);
            num_variables = size(x_diff,1);
            x_diff_weighted = zeros(num_variables, num_sigma_points-1);
            for iPoints = 2:num_sigma_points
                x_diff_weighted(:,iPoints-1) = ...
                    sqrt(weights_cov(1,iPoints))*x_diff(:,iPoints);
            end
            M = horzcat(x_diff_weighted, square_root_covariance);
            [foo, S] = qr(M.');
            S_trimmed = S(1:num_variables,1:num_variables);
            if sign(weights_cov(1,1)) == 1
                operation = '+';
            elseif sign(weights_cov(1,1)) == -1
                operation = '-';
            else
                error('QR decomposition error');
            end
            % output = cholupdate(S_trimmed, sqrt(norm(weights_cov(1,1)))*x_diff(:,1), operation);
            hoge = cholupdate(S_trimmed, sqrt(norm(weights_cov(1,1)))*x_diff(:,1), operation);
            output = hoge.';
        end

        function output = qr_decomposition_pure(this, A)
            [Q, R] = qr(A.');
            column_size = size(R,2);
            R_tilde = R(1:column_size, 1:column_size);
            hoge = R_tilde.';
            output = hoge;
        end

    end
end