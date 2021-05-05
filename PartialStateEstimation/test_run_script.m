clc
clear
close all
addpath(genpath('../../state_estimation'));
import matlab.unittest.TestSuite
suiteClass = TestSuite.fromClass(?TestPartialStateEstimation);
result = run(suiteClass);