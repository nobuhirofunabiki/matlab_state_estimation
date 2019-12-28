classdef PF_LinearDynamics < ParticleFilterBase
    properties (SetAccess = private)
        discrete_system_matrix % Discrete matrix based on the same time step with the estimation period
        process_noise_covmat
        direct_roughening_covmat
        num_dimensions
    end

    methods (Access = protected)
        function obj = PF_LinearDynamics(args)
            obj@ParticleFilterBase(args);
            obj.checkConstructorArguments(args);
            obj.num_dimensions = args.num_dimensions;
            obj.discrete_system_matrix = args.discrete_system_matrix;
            obj.process_noise_covmat = args.process_noise_covmat;
            obj.direct_roughening_covmat = args.direct_roughening_covmat;
        end
    end

    methods (Access = private)
        function checkConstructorArguments(this, args)
            disp("Check constructor arguments for PF_LinearDynamics");
            num_vars = this.number_variables;
            assert(isequal(size(args.discrete_system_matrix), [num_vars, num_vars]), ...
                "discrete_system_matrix is NOT correct size matrix");
            assert(isequal(size(args.process_noise_covmat), [num_vars, num_vars]), ...
                "process_noise_covmat is NOT correct size matrix");
            assert(isequal(size(args.direct_roughening_covmat), [num_vars, num_vars]), ...
                "direct_roughening_covmat is NOT correct size matrix");
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

    methods (Access = public)
        % Setters -----------------------------------------
        function output = getNumberDimensions(this)
            output = this.num_dimensions;
        end
    end

end