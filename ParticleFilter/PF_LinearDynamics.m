classdef PF_LinearDynamics < ParticleFilterBase
    properties (SetAccess = private)
        discrete_system_matrix % Discrete matrix based on the same time step with the estimation period
        process_noise_covmat
        direct_roughening_covmat
        num_dimensions
    end

    methods
        function obj = PF_LinearDynamics(args)
            obj@ParticleFilterBase(args);
            obj.discrete_system_matrix = args.discrete_system_matrix;
            obj.process_noise_covmat = args.process_noise_covmat;
            obj.direct_roughening_covmat = args.direct_roughening_covmat;
            obj.num_dimensions = args.num_dimensions;
        end

        function output = propagateParticleStates(this)
            for iParticles = 1:this.number_particles
                this.particle_states(:,iParticles) = ...
                    this.discrete_system_matrix * this.prior_particle_states(:,iParticles) ...
                    + sqrt(this.process_noise_covmat) * randn(length(this.particle_states(:,iParticles)),1) ...
                    + sqrt(this.direct_roughening_covmat) * randn(length(this.particle_states(:,iParticles)),1);
            end
        end

        % Setters -----------------------------------------
        function output = getNumberDimensions(this)
            output = this.num_dimensions;
        end

    end
end