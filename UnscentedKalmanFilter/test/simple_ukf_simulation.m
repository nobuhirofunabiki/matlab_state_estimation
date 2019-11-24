clear
clc
close all

restoredefaultpath
addpath(genpath('../../UnscentedKalmanFilter'));

rng('default');

% States
position_true = [0; 0];
position_est = [1; 1];
position_cov = [10 0; 0 10];

% Dynamics
A = [1 0; 0 1];
sys_cov = [0.1 0; 0 0.1];

% Measurements
range_list = [...
5   5   -5  -5;
-5  5   5   -5];
range_cov_diag = 0.3*ones(4,1);
range_cov = diag(range_cov_diag);

% Instanciation of Unscented Kalman Filter
% args_ukf.alpha = 1e-3;
args_ukf.alpha = 0.5;
args_ukf.beta = 2.0;
args_ukf.kappa = 0;
args_ukf.state_dim = 2;
args_ukf.state_est = position_est;
args_ukf.state_cov = position_cov;
args_ukf.range_list = range_list;
args_ukf.obs_cov = range_cov;
args_ukf.A = A;
args_ukf.sys_cov = sys_cov;
ukf_ = UkfTargetTrackingRange(args_ukf);

% Simulation
for iSteps = 1:100
    
end


ukf_.executeUnscentedKalmanFilter();