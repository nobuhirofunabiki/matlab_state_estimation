classdef TestPartialStateEstimation < matlab.unittest.TestCase
    properties
        object_
    end

    methods(TestMethodSetup)
        function instanciatePartialStateEstimation(testCase)
            args.num_agents = 5;
            args.num_dimensions = 1;
            args.num_variables = 2 * args.num_agents * args.num_dimensions;
            args.agent_id = 1;
            args.process_noise_covmat = zeros(args.num_variables, args.num_variables);
            args.state_vector = zeros(args.num_variables, 1);
            args.sigma_position = 0;
            args.sigma_velocity = 0;
            testCase.object_ = PartialStateEstimation(args);
        end
    end

    methods(Test)
        function testCheckNetworkUpdate(testCase)
            object_ = testCase.object_;
            old_adjacent_matrix = [...
                0 1 1 0 0;
                1 0 0 0 0;
                1 0 0 1 0;
                0 0 1 0 1;
                0 0 0 1 0];
            new_adjacent_matrix = [...
                0 1 0 0 0;
                1 0 1 0 0;
                0 1 0 0 1;
                0 0 0 0 1;
                0 0 1 1 0];
            agent_id = object_.getAgentID();
            estimated_agents = old_adjacent_matrix(agent_id,:);
            estimated_agents(1,agent_id) = 1;
            object_.setCurrentEstimatedAgents(estimated_agents);
            object_.checkNetworkUpdate(new_adjacent_matrix);
            result_previous_estimated_agents = object_.getPreviousEstimatedAgents();
            result_current_estimated_agents = object_.getCurrentEstimatedAgents();
            result_number_estimated = object_.getNumberEstimated();
            expected_previous_estimated_agents = [1 1 1 0 0];
            expected_current_estimated_agents = [1 1 0 0 0];
            expected_number_estimated = 2;
            testCase.verifyEqual(result_previous_estimated_agents, expected_previous_estimated_agents);
            testCase.verifyEqual(result_current_estimated_agents, expected_current_estimated_agents);
            testCase.verifyEqual(result_number_estimated, expected_number_estimated);
        end

        function testMakeStateCovarianceMatrix1(testCase)
            object_ = testCase.object_;
            shared_state_covmat = zeros(10, 10);
            for iAgents = 1:5
                for jAgents = 1:5
                    shared_state_covmat(2*(iAgents-1)+1:2*iAgents, 2*(jAgents-1)+1:2*jAgents) ...
                        = (10*iAgents+jAgents)*ones(2,2);
                end
            end
            previous_state_covmat = zeros(8,8);
            for iAgents = 1:4
                for jAgents = 1:4
                    previous_state_covmat(2*(iAgents-1)+1:2*iAgents, 2*(jAgents-1)+1:2*jAgents) ...
                        = (10*iAgents+jAgents)*ones(2,2);
                end
            end
            previous_estimated_agents = [1 1 1 1 0];
            current_estimated_agents = [1 1 1 0 1];
            num_estimated = sum(current_estimated_agents);
            object_.setPreviousEstimatedAgents(previous_estimated_agents);
            object_.setCurrentEstimatedAgents(current_estimated_agents);
            object_.setNumberEstimated(num_estimated);
            object_.setStateCovarianceMatrix(previous_state_covmat);
            result = object_.makeStateCovarianceMatrix(shared_state_covmat);
            expected = zeros(8,8);
            for iAgents = 1:3
                for jAgents = 1:3
                    expected(2*(iAgents-1)+1:2*iAgents, 2*(jAgents-1)+1:2*jAgents) ...
                        = (10*iAgents+jAgents)*ones(2,2);
                end
            end
            expected(7:8,7:8) = 55*ones(2,2);
            testCase.verifyEqual(result, expected);
        end

        function testMakeStateCovarianceMatrix2(testCase)
            object_ = testCase.object_;
            shared_state_covmat = zeros(10, 10);
            for iAgents = 1:5
                for jAgents = 1:5
                    shared_state_covmat(2*(iAgents-1)+1:2*iAgents, 2*(jAgents-1)+1:2*jAgents) ...
                        = (10*iAgents+jAgents)*ones(2,2);
                end
            end
            previous_state_covmat = zeros(6,6);
            for iAgents = 1:3
                for jAgents = 1:3
                    previous_state_covmat(2*(iAgents-1)+1:2*iAgents, 2*(jAgents-1)+1:2*jAgents) ...
                        = (10*iAgents+jAgents)*ones(2,2);
                end
            end
            previous_estimated_agents = [1 1 1 0 0];
            current_estimated_agents = [0 0 0 1 1];
            num_estimated = sum(current_estimated_agents);
            object_.setPreviousEstimatedAgents(previous_estimated_agents);
            object_.setCurrentEstimatedAgents(current_estimated_agents);
            object_.setNumberEstimated(num_estimated);
            object_.setStateCovarianceMatrix(previous_state_covmat);
            result = object_.makeStateCovarianceMatrix(shared_state_covmat);
            expected = zeros(4,4);
            expected(1:2,1:2) = 44*ones(2,2);
            expected(3:4,3:4) = 55*ones(2,2);
            testCase.verifyEqual(result, expected);
        end
    end
end