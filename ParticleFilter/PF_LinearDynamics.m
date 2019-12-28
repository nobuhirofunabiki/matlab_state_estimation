classdef PF_LinearDynamics < ParticleFilterBase

    properties (Abstract = true, SetAccess = protected)
        discrete_system_matrix
        process_noise_covmat
        direct_roughening_covmat
    end

    methods (Access = protected)
        function obj = PF_LinearDynamics(args)
            obj@ParticleFilterBase(args);
        end
    end

    methods (Access = protected)
        function output = propagateParticleStates(this)
            for iParticles = 1:this.number_particles
                this.particle_states(:,iParticles) = ...
                    this.discrete_system_matrix * this.prior_particle_states(:,iParticles) ...
                    + sqrt(this.process_noise_covmat) * randn(length(this.particle_states(:,iParticles)),1) ...
                    + sqrt(this.direct_roughening_covmat) * randn(length(this.particle_states(:,iParticles)),1);
            end
        end
    end

end