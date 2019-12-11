classdef ParticleFilterTargetTrackingRange < ParticleFilterRangePosition
    properties (SetAccess = protected)
        range_list
    end

    methods
        function obj = ParticleFilterTargetTrackingRange(args)
            obj@ParticleFilterRangePosition(args);
            obj.range_list = args.range_list;
        end

        function output = calculateEstimatedMeasurements(this, args)
            iParticles = args.iParticles;
            NUM_DIMS = args.number_dimensions;
            NUM_MEASURES = size(this.range_list,2);
            est_measures = zeros(NUM_MEASURES,1);
            position = this.particle_states(1:NUM_DIMS, iParticles);
            for iMeasures = 1:NUM_MEASURES
                est_measures(iMeasures,1) = ...
                    norm(position - this.range_list(:,iMeasures));
            end
            output = est_measures;
        end
    end

end