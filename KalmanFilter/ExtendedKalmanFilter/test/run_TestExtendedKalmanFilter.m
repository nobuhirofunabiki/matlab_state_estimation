addpath(genpath('../../../state_estimation'));
import matlab.unittest.TestSuite
suiteClass = TestSuite.fromClass(?TestExtendedKalmanFilter);
result = run(suiteClass);