classdef PF_TargetTrackingRangeLandmarks < PF_LinearDynamics
    

    properties (SetAccess = protected)
        range_sensor_       % Class instance of RangeMeasurementLandmarks
        counter
        discrete_system_matrix
        process_noise_covmat
        direct_roughening_covmat
    end
    properties (SetAccess = immutable)
        num_dimensions
        num_agents
        init_state_vector
        init_state_covmat
    end

    methods (Access = public)
        function obj = PF_TargetTrackingRangeLandmarks(args)
            obj@PF_LinearDynamics(args);
            obj.checkConstructorArguments(args);
            obj.range_sensor_               = RangeMeasurementLandmarks(args.rml);
            obj.counter                     = 0;
            obj.discrete_system_matrix      = args.discrete_system_matrix;
            obj.process_noise_covmat        = args.process_noise_covmat;
            obj.direct_roughening_covmat    = args.direct_roughening_covmat;
            obj.num_dimensions              = args.num_dimensions;
            obj.num_agents                  = args.num_agents;
            obj.init_state_vector           = args.init_state_vector;
            obj.init_state_covmat           = args.init_state_covmat;
            obj.setInitialParticleStates();
        end
    end

    methods (Access = private)
        function checkConstructorArguments(this, args)
            disp("Check constructor arguments for PF_TargetTrackingRangeLandmarks");
            assert(isequal(size(args.init_state_covmat), [this.number_variables, this.number_variables]), ...
                "init_state_covmat is NOT correct size matrix");
        end
    end

    methods (Access = public)
        function executeParticleFiltering(this, measurements)
            this.updateParticles(measurements);
            this.resampleParticles();
            % TODO: This resampling for velocity should be implemented in a more logical way
            if (this.counter < 10)
                this.setParticleStatesOnlyVelocity();
            end
            this.counter = this.counter + 1;
            this.computeEstimatedStates();
            this.prepareForNextFiltering();
        end
    end

    methods (Access = protected)
        function updateParticleWeights(this, measurements)
            for iParticles = 1:this.number_particles
                NUM_DIMS = this.num_dimensions;
                position = this.particle_states(1:NUM_DIMS, iParticles);
                this.range_sensor_.computeMeasurementVector(position, false);
                predicted_measurements = this.range_sensor_.getMeasurements();
                diff_measurements = transpose(measurements - predicted_measurements);
                obs_covmat = this.range_sensor_.getMeasureCovarinaceMatrix();
                for iMeasures = 1:length(diff_measurements)
                    noise_sigma = sqrt(obs_covmat(iMeasures, iMeasures));
                    likelihood = normpdf(diff_measurements(1,iMeasures), 0, noise_sigma);
                    this.weights(iParticles,1) = this.weights(iParticles,1) * likelihood;
                end
            end
        end

        function setParticleStatesOnlyVelocity(this)
            num_dims = this.num_dimensions;
            for iParticles = 1:this.number_particles
                particle_states = mvnrnd(this.particle_states(:,iParticles), this.init_state_covmat);
                for iAgents = 1:this.num_agents
                    for iDims = 1:num_dims
                        this.particle_states(2*num_dims*(iAgents-1)+num_dims+iDims, iParticles) ...
                            = particle_states(2*num_dims*(iAgents-1)+num_dims+iDims);
                    end
                end
            end
        end

    end
end