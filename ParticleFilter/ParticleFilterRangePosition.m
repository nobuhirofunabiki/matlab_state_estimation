classdef ParticleFilterRangePosition < ParticleFilterBase
    properties (SetAccess = protected)
        system_matrix % Discrete matrix based on the same time step with the estimation period
        sys_covmat
        obs_covmat
        num_measures
        num_agents
        num_dims
    end

    methods
        function obj = ParticleFilterRangePosition(args)
            obj@ParticleFilterBase(args);
            obj.num_measures = args.num_measures;
            obj.num_dims = args.num_dims;
            obj.num_agents = args.num_agents;
            obj.system_matrix = zeros(obj.number_variables, obj.number_variables);
            obj.sys_covmat = zeros(size(obj.system_matrix));
            covmat_diag = 10.^10*ones(obj.num_measures,1);
            obj.obs_covmat = diag(covmat_diag);
        end

        function executeParticleFiltering(this, args)
            args_updateParticles.time_step = args.time_step;
            args_updateParticles.measurements = args.measurements;
            args_updateParticles.number_agents = args.number_agents;
            args_updateParticles.number_dimensions = args.number_dimensions;
            this.updateParticles(args_updateParticles);
            this.resampleParticles();
            this.computeExtimatedStates();
            this.prepareForNextFiltering();
        end

        % Setters -------------------------------------------------

        function setFirstParticleStates(this, initial_estimates)
            for iParticles = 1:this.number_particles
                this.prior_particle_states(:,iParticles) = ...
                    mvnrnd(initial_estimates, this.sys_covmat);
            end
        end

        function setObservationCovarianceMatrix(this, covmat)
            if (size(covmat) ~= size(this.obs_covmat))
                error("Invalid size of measurement covariance matrix")
            end
            this.obs_covmat = covmat;
        end

        function setStateCovarianceMatrix(this, args)
            P = zeros(size(this.sys_covmat));
            num_dims = this.num_dims;
            for iAgents = 1:this.num_agents
                for iDims = 1:num_dims
                    P(2*num_dims*(iAgents-1)+iDims, 2*num_dims*(iAgents-1)+iDims) = args.position_sigma^2;
                    P(2*num_dims*(iAgents-1)+num_dims+iDims, 2*num_dims*(iAgents-1)+num_dims+iDims) = args.velocity_sigma^2;
                end
            end
            this.sys_covmat = P;
        end

        function setDiscreteSystemMatrix(this, discrete_system_matrix)
            num_agents = this.num_agents;
            num_vars = this.number_variables;
            num_dims = this.num_dims;
            Ad = zeros(num_vars, num_vars);
            for iAgents = 1:num_agents
                Ad(1+2*num_dims*(iAgents-1):2*num_dims*iAgents,...
                    1+2*num_dims*(iAgents-1):2*num_dims*iAgents)...
                = discrete_system_matrix;
            end
            this.system_matrix = Ad;
        end
    end

    methods
        function updateParticles(this, args)
            args_temp.number_agents = args.number_agents;
            args_temp.number_dimensions = args.number_dimensions;
            for iParticles = 1:this.number_particles
                % Update the particle states
                % this.particle_states(:,iParticles) = ...
                %     this.eq_sys(args.time_step, ...
                %         this.prior_particle_states(:,iParticles), ...
                %         this.generator_sys_noise());
                this.propagateParticleState(iParticles);
                % Update the importance factors
                args_temp.iParticles = iParticles;
                estimated_measurements = calculateEstimatedMeasurements(args_temp);
                means = zeros(this.num_measures,1);
                likelihood = mvnpdf(args.measurements - estimated_measurements, means, this.obs_covmat);
                % this.weights(iParticles,1) = ...
                %     this.prior_weights(iParticles,1) ...
                %     + this.pdf_likelihood(args.time_step, ...
                %         args.measurements, estimated_measurements);
                this.weights(iParticles,1) = this.prior_weights(iParticles,1) + likelihood;
            end
            % Normalization of the importance factors
            this.weights = this.weights./sum(this.weights);
        end

        function output = propagateParticleState(this, iParticles)
            % Only for the linear system equation
            this.particle_states(:,iParticles) = ...
                this.system_matrix * this.prior_particle_states(:,iParticles) ...
                + mvnrnd(zeros(size(this.particle_states(:,iParticles))), this.sys_covmat);
        end

        function output = calculateEstimatedMeasurements(this, args)
            NUM_AGENTS = args.number_agents;
            NUM_DIMS = args.number_dimensions;
            NUM_MEASURES = nchoosek(NUM_AGENTS,2) + NUM_DIMS*NUM_AGENTS;
            iParticles = args.iParticles;
            est_measures = zeros(NUM_MEASURES, 1);
            iMeasures = 0;
            % Estimated range measurements
            for iAgents = 1:NUM_AGENTS-1
                for jAgents = iAgents+1:NUM_AGENTS
                    iMeasures = iMeasures + 1;
                    est_pos_iAgent = this.particle_states(2*NUM_DIMS*(iAgents-1)+1:2*NUM_DIMS*(iAgents-1)+2, iParticles);
                    est_pos_jAgent = this.particle_states(2*NUM_DIMS*(jAgents-1)+1:2*NUM_DIMS*(jAgents-1)+2, iParticles);
                    est_measures(iMeasures, 1) = norm(est_pos_iAgent - est_pos_jAgent);
                end
            end
            % Estimated position measurements
            for iAgents = 1:NUM_AGENTS
                for iDims = 1:NUM_DIMS
                    iMeasures = iMeasures + 1;
                    est_mesaures(iMeasures) = this.particle_states(2*NUM_DIMS*(iAgents-1)+iDims);
                end
            end
            output = est_measures;
        end
    end
end