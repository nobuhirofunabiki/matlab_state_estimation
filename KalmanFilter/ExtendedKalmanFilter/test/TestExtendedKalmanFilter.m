classdef TestExtendedKalmanFilter < matlab.unittest.TestCase
    properties
        ekf_
    end

    methods (TestMethodSetup)
        function createExtendedKalmanFilterInstance(testCase)
            NUM_VARIABLE = 6;
            NUM_MEASURE = 12;
            args.num_variable = NUM_VARIABLE;
            args.num_measure = NUM_MEASURE;
            testCase.ekf_ = ExtendedKalmanFilter(args);
            covariance_matrix = [...
                1 2 3 4 5 6;
                7 8 9 10 11 12;
                13 14 15 16 17 18;
                19 20 21 22 23 24;
                25 26 27 28 29 30;
                31 32 33 34 35 36];
            testCase.ekf_.setStateCovarianceMatrix(covariance_matrix);
        end
    end

    methods (Test)
        function testGetPartialStateCovarianceMatrix(testCase)
            ekf_ = testCase.ekf_;
            partial_covmat_actual = ...
                ekf_.getPartialStateCovarianceMatrix(2,4,3,5);
            partial_covmat_exptected = [...
                9 10 11;
                15 16 17;
                21 22 23];
            testCase.verifyEqual(partial_covmat_actual, partial_covmat_exptected);
        end
    end
end