classdef ParticleFilterBase < handle
    properties (SetAccess = protected)
        number_particles
        number_variables
        estimated_state;
        weights
        particle_states
        prior_particle_states
        resample_percentage
    end
    properties (Abstract = true, SetAccess = immutable)
        init_state_vector
        init_state_covmat
    end

    methods (Access = protected)
        function obj = ParticleFilterBase(args)
            obj.number_particles = args.number_particles;
            obj.number_variables = args.number_variables;
            obj.estimated_state = zeros(obj.number_variables, 1);
            obj.weights = zeros(args.number_particles, 1);
            for iParticles = 1:obj.number_particles
                obj.weights(iParticles, 1) = 1./obj.number_particles;
            end
            obj.particle_states = ...
                zeros(args.number_variables, args.number_particles);
            obj.prior_particle_states = zeros(size(obj.particle_states));
            obj.resample_percentage = args.resample_percentage;
        end
    end

    methods (Abstract = true, Access = protected)
        propagateParticleStates(this);
        updateParticleWeights(this);
    end

    methods (Access = public)
        function executeParticleFiltering(this, measurements)
            this.updateParticles(measurements);
            this.resampleParticles();
            this.computeEstimatedStates();
            this.prepareForNextFiltering();
        end

        % Getters -------------------------------------------------------
        function output = getStateVector(this)
            output = this.estimated_state;
        end

        function output = getParticleNumber(this);
            output = this.number_particles;
        end

        function output = getParticleStates(this)
            output = this.particle_states;
        end
    end

    methods (Access = protected)
        function updateParticles(this, measurements)
            this.propagateParticleStates();
            this.updateParticleWeights(measurements);
            this.normarizeWeights();
        end

        function resampleParticles(this)
            Ns = this.number_particles;
            effective_sample_size = 1/sum(this.weights.^2);
            sample_threshold = Ns * this.resample_percentage;
            if effective_sample_size < sample_threshold
                edges = min([0 cumsum(this.weights)'],1); % protect against accumulated round-off
                edges(end) = 1;                 % get the upper edge exact
                u1 = rand/Ns;
                % this works like the inverse of the empirical distribution and returns
                % the interval where the sample is to be found
                [~, idx] = histc(u1:1/Ns:1, edges);
                % idx = randsample(1:Ns, Ns, true, this.weights);
                particle_states = this.particle_states;
                for iParticles = 1:Ns
                    this.particle_states(:,iParticles) = particle_states(:,idx(1,iParticles));
                end
                this.weights = transpose(repmat(1/Ns, 1, Ns));
            end
        end

        function normarizeWeights(this)
            sum_weghts = sum(this.weights);
            if sum_weghts ~= 0
                this.weights = this.weights./sum(this.weights);
            else
                this.weights = zeros(this.number_particles, 1) + 1./this.number_particles;
            end
        end

        function computeEstimatedStates(this)
            updated_states = zeros(this.number_variables, 1);
            for iParticles = 1:this.number_particles
                updated_states = updated_states ...
                    + this.weights(iParticles,1) * this.particle_states(:,iParticles);
            end
            this.estimated_state = updated_states;
        end

        function prepareForNextFiltering(this)
            this.prior_particle_states = this.particle_states;
        end

        % Setters -------------------------------------------------------

        function setInitialParticleStates(this)
            for iParticles = 1:this.number_particles
                this.prior_particle_states(:,iParticles) = ...
                    mvnrnd(this.init_state_vector, this.init_state_covmat);
            end
        end

    end
end