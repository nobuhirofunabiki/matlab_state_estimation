%% clear memory, screen, and close all figures
clear
clc
% close all;

restoredefaultpath
addpath(genpath('../../ParticleFilter'));

rng('default');

%% Setup ----------------------------------------------------------

% Definition of system and observation equations
eq_sys = @(k, prior_state, state_noise) ...
    prior_state/2 + 25*prior_state/(1+prior_state^2) + 8*cos(1.2*k) + state_noise;
eq_obs = @(k, state, obs_noise) ...
    state^2/20 + obs_noise;

% PDF of system noise and noise generator function
sigma_sys_noise = sqrt(10);
pdf_sys_noise = @(u) normpdf(u, 0, sigma_sys_noise);
generator_sys_noise = @(u) normrnd(0, sigma_sys_noise);

% PDF of observation noise and noise gerator function
sigma_obs_noise = sqrt(1);
pdf_obs_noise = @(u) normpdf(u, 0, sigma_obs_noise);
generator_obs_noise = @(u) normrnd(0, sigma_obs_noise);

% Initial PDF
gen_x0 = @(x) normrnd(0, sqrt(10));

% Observation likelihood PDF p(y[k] | x[k])
% (under the suposition of additive process noise)
pdf_likelihood = @(k, yk, xk) pdf_obs_noise(yk - eq_obs(k, xk, 0));

% Instanciation of particle filter
args_pf.number_particles = 200;
args_pf.number_variables = 1;
args_pf.prior_particle_state = gen_x0;
args_pf.resample_percentage = 0.50;
args_pf.eq_sys = eq_sys;
args_pf.eq_obs = eq_obs;
args_pf.pdf_likelihood = pdf_likelihood;
args_pf.generator_sys_noise = generator_sys_noise;
particle_filter_ = ParticleFilter(args_pf);

TIME_STEP_TOTAL = 40;
time_list = zeros(1,TIME_STEP_TOTAL);

% Visualization
% Visualization: True system
viz_state = zeros(1,TIME_STEP_TOTAL);
viz_measures = zeros(1,TIME_STEP_TOTAL);
viz_sys_noise = zeros(1,TIME_STEP_TOTAL);
viz_obs_noise = zeros(1,TIME_STEP_TOTAL);
% Visualization: Estimated system
viz_state_est = zeros(1,TIME_STEP_TOTAL);

%% Simulation -------------------------------------------------------

% Set initial values
time_list(:,1) = 1;
sys_noise = 0;
obs_noise = generator_obs_noise(sigma_obs_noise);
state = 0;
measurement = eq_obs(1, state, obs_noise);

viz_state(:,1) = state;
viz_measures(:,1) = measurement;
viz_sys_noise(:,1) = sys_noise;
viz_obs_noise(:,1) = obs_noise;

% Iterations
for iCounts = 2:TIME_STEP_TOTAL
    time_list(:,iCounts) = iCounts;
    sys_noise = generator_sys_noise();
    obs_noise = generator_obs_noise();
    state = eq_sys(iCounts, state, sys_noise);
    measurement = eq_obs(iCounts, state, obs_noise);

    args_exe_pf.time_step = iCounts;
    args_exe_pf.measurements = measurement;
    particle_filter_.executeParticleFiltering(args_exe_pf);

    viz_state(:,iCounts) = state;
    viz_measures(:,iCounts) = measurement;
    viz_sys_noise(:,iCounts) = sys_noise;
    viz_obs_noise(:,iCounts) = obs_noise;
    viz_state_est(:,iCounts) = particle_filter_.getStateVector();
end

%% Visualization

figure
plot(time_list, viz_state, 'k', time_list, viz_state_est, 'b');
legend('true state', 'estimated_state')
grid on
