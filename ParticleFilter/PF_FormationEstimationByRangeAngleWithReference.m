classdef PF_FormationEstimationByRangeAngleWithReference < ...
    PF_LinearDynamics & ...
    MultiagentUtilityBase

    properties (SetAccess = protected)
        discrete_system_matrix
        process_noise_covmat
        direct_roughening_covmat
    end
    properties (SetAccess = private)
        range_sensor_       % Class instance of RangeMeasurementMultiAgentWithReference
        angle_sensor_       % Class instance of AngleMeasurementMultiAgentWithReference
        counter
    end
    properties (SetAccess = immutable)
        num_dimensions
        num_agents
        init_state_vector
        init_state_covmat
    end

    methods (Access = public)
        function obj = PF_FormationEstimationByRangeAngleWithReference(args)
            obj@PF_LinearDynamics(args);
            obj.checkConstructorArguments(args);
            obj.num_dimensions              = args.num_dimensions;
            obj.num_agents                  = args.num_agents;
            obj.range_sensor_               = RangeMeasurementMultiAgentWithReference(args.range_sensor);
            obj.angle_sensor_               = AngleMeasurementMultiAgentWithReference(args.angle_sensor);
            obj.counter                     = 0;
            obj.discrete_system_matrix      = obj.createDiscreteSystemMatrix(args.discrete_system_matrix);
            obj.process_noise_covmat        = args.process_noise_covmat;
            obj.direct_roughening_covmat    = args.direct_roughening_covmat;
            init_state                      = obj.createInitialStateVectorAndCovariance(args);
            obj.init_state_vector           = init_state.init_state_vector;
            obj.init_state_covmat           = init_state.init_state_covmat;
            obj.setInitialParticleStates();
        end
    end

    methods (Access = private)
        function checkConstructorArguments(this, args)
            disp("Check constructor arguments for PF_LinearDynamics");
            num_vars = this.number_variables;
            num_dims = args.num_dimensions;
            assert(isequal(size(args.discrete_system_matrix), [2*num_dims, 2*num_dims]), ...
                "discrete_system_matrix is NOT a correct size matrix");
            assert(isequal(size(args.process_noise_covmat), [num_vars, num_vars]), ...
                "process_noise_covmat is NOT a correct size matrix");
            assert(isequal(size(args.direct_roughening_covmat), [num_vars, num_vars]), ...
                "direct_roughening_covmat is NOT a correct size matrix");
            assert(isequal(size(args.init_state_vector), [num_vars, 1]), ...
                "init_state_vector is NOT a correct size vector");
        end
    end

    methods (Access = public)
        function executeParticleFiltering(this, measurements, adjacent_matrix, position_ref)
            this.updateParticles(measurements, adjacent_matrix, position_ref);
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
        function updateParticles(this, measurements, adjacent_matrix, position_ref)
            this.propagateParticleStates();
            this.updateParticleWeights(measurements, adjacent_matrix, position_ref);
            this.normarizeWeights();
        end

        function updateParticleWeights(this, measurements, adjacent_matrix, position_ref)
            % adjacent_matrix should have 'range' field and 'angle' field
            for iParticles = 1:this.number_particles
                NUM_DIMS = this.num_dimensions();
                positions = this.getPositionVector(iParticles);

                % Range measurements
                this.range_sensor_.computeMeasurementVector(positions, position_ref, false);
                this.range_sensor_.updateMeasurementCovarianceMatrix(adjacent_matrix.range);
                predicted_range_measurements = this.range_sensor_.getMeasurements();
                diff_range_measurements = transpose(measurements.ranges - predicted_range_measurements);
                obs_covmat_range = this.range_sensor_.getMeasureCovarinaceMatrix();
                for iMeasures = 1:length(diff_range_measurements)
                    noise_sigma = sqrt(obs_covmat_range(iMeasures, iMeasures));
                    likelihood = normpdf(diff_range_measurements(1,iMeasures), 0, noise_sigma);
                    this.weights(iParticles,1) = this.weights(iParticles,1) * likelihood;
                end

                % Angle measurements
                this.angle_sensor_.computeMeasurementVector(positions, position_ref, false);
                this.angle_sensor_.setMeasurementCovarianceMatrix(adjacent_matrix.angle);
                predicted_angle_measurements = this.angle_sensor_.getMeasurements();
                % TODO: How to tuckle the following singular point problem?
                diff_measures = predicted_angle_measurements - measurements.angles;
                for iMeasures = 1:length(diff_measures)
                    if (abs(diff_measures(iMeasures,1)) >= pi)
                        if (predicted_angle_measurements(iMeasures,1) > measurements.angles(iMeasures,1))
                            predicted_angle_measurements(iMeasures,1) = predicted_angle_measurements(iMeasures,1) - 2*pi;
                        else
                            predicted_angle_measurements(iMeasures,1) = predicted_angle_measurements(iMeasures,1) + 2*pi;
                        end
                    end
                end
                diff_angle_measurements = transpose(measurements.angles - predicted_angle_measurements);
                obs_covmat_angle = this.angle_sensor_.getMeasureCovarinaceMatrix();
                for iMeasures = 1:length(diff_angle_measurements)
                    noise_sigma = sqrt(obs_covmat_angle(iMeasures, iMeasures));
                    likelihood = normpdf(diff_angle_measurements(1,iMeasures), 0, noise_sigma);
                    this.weights(iParticles,1) = this.weights(iParticles,1) * likelihood;
                end

            end
        end

        function output = createInitialStateVectorAndCovariance(this, args)
            sigma_position = args.init_sigma_position;
            sigma_velocity = args.init_sigma_velocity;
            output.init_state_vector = args.init_state_vector;
            output.init_state_covmat = ...
                this.createStateCovarianceMatrix(sigma_position, sigma_velocity);
        end
    end

    methods (Access = private)
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

        function output = getPositionVector(this, iParticles)
            num_dims = this.num_dimensions;
            output = zeros(num_dims*this.num_agents,1);
            for iAgents = 1:this.num_agents
                output(num_dims*(iAgents-1)+1:num_dims*iAgents,1) ...
                    = this.particle_states(2*num_dims*(iAgents-1)+1:2*num_dims*(iAgents-1)+num_dims,iParticles);
            end
        end
    end
    
end