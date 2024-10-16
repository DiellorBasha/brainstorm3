classdef graph_unit_tests < matlab.unittest.TestCase
% Copyright 2021 The MathWorks, Inc.
    methods (Test)
        function check_invalid_start_1(testCase)
            adjMatrix = graph_unit_tests.graph_straight_seq();
            startIdx = -1;
            endIdx = 2;
            expOut = -199;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Invalid start index, idx<1');
        end

        function check_invalid_start_2(testCase)
            adjMatrix = graph_unit_tests.graph_straight_seq();
            startIdx = 12;
            endIdx = 2;
            expOut = -99;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Invalid start index, idx>NodeCnt');
        end

        function check_invalid_end_1(testCase)
            adjMatrix = graph_unit_tests.graph_straight_seq();
            startIdx = 1;
            endIdx = -3;
            expOut = -199;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Invalid end index, idx<1');
        end

        function check_invalid_end_2(testCase)
            adjMatrix = graph_unit_tests.graph_straight_seq();
            startIdx = 1;
            endIdx = 12;
            expOut = -99;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Invalid end index, idx>NodeCnt');
        end

        function check_longest_path(testCase)
            adjMatrix = graph_unit_tests.graph_straight_seq();
            startIdx = 1;
            endIdx = 4;
            expOut = 3;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Longest theoretic path');
        end

        function check_unity_path(testCase)
            adjMatrix = graph_unit_tests.graph_all_edge();
            startIdx = 2;
            endIdx = 3;
            expOut = 1;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Path length 1');
        end

        function check_non_unique(testCase)
            adjMatrix = graph_unit_tests.graph_square();
            startIdx = 4;
            endIdx = 2;
            expOut = 2;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Non-unique path');
        end

        function check_no_path(testCase)
            adjMatrix = graph_unit_tests.graph_disconnected_components();
            startIdx = 1;
            endIdx = 5;
            expOut = -1;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'No path');
        end

        function check_edgeless_graph(testCase)
            adjMatrix = zeros(20,20);
            startIdx = 1;
            endIdx = 18;
            expOut = -1;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Edgeless graph');
        end

        function check_edgeless_start(testCase)
            adjMatrix = graph_unit_tests.graph_some_nodes_edgeless();
            startIdx = 1;
            endIdx = 4;
            expOut = -1;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Edgeless graph');
        end

        function check_edgeless_end(testCase)
            adjMatrix = graph_unit_tests.graph_some_nodes_edgeless();
            startIdx = 3;
            endIdx = 1;
            expOut = -1;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Edgeless graph');
        end

        
        function check_edgeless_graph_self_loop(testCase)
            adjMatrix = zeros(20,20);
            startIdx = 16;
            endIdx = 16;
            expOut = 0;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Self loop in edgeless graph');
        end

        function check_start_end_same(testCase)
            adjMatrix = graph_unit_tests.graph_all_edge();
            startIdx = 3;
            endIdx = 3;
            expOut = 0;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Start and end index same');
        end

        function check_invalid_idx_empty_adj(testCase)
            adjMatrix = [];
            startIdx = 1;
            endIdx = 1;
            expOut = -99;
            verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, 'Degenerate empty graph');
        end

        function check_invalid_nonsquare(testCase)
    adjMatrix = zeros(2,3);
    startIdx = 1;
    endIdx = 1;
    expOut = -9;
    verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, ...
        'Graph is not square');
end

function check_invalid_entry(testCase)
    adjMatrix = 2*ones(4,4);
    startIdx = 1;
    endIdx = 1;
    expOut = -9;
    verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, ...
        'Adjacency matrix is not valid');
end

function check_invalid_noninteger_startnode(testCase)
    adjMatrix = zeros(4,4);
    startIdx = 1.2;
    endIdx = 1;
    expOut = -19;
    verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, ...
        'Start node is not an integer');
end

function check_invalid_noninteger_endnode(testCase)
    adjMatrix = zeros(4,4);
    startIdx = 1;
    endIdx = 2.2;
    expOut = -29;
    verify_path_length(testCase, adjMatrix, startIdx, endIdx, expOut, ...
        'End node is not an integer');
end
    end

    methods (Static)
        % Utility functions to create common adjacency graph matrices
        function adj = graph_straight_seq()
            % Create the graph:
            % 1---2---3---4

            adj = [0 1 0 0; ...
                1 0 1 0; ...
                0 1 0 1; ...
                0 0 1 0];
        end

        function adj = graph_square()
            % Create the graph:
            %   1---2
            %   |   |
            %   4---3        

            adj = [0 1 0 1; ...
                1 0 1 0; ...
                0 1 0 1; ...
                1 0 1 0];
        end

        function adj = graph_all_edge()
            % Create the graph:
            %   1---2
            %   |\ /|
            %   |/ \|
            %   4---3

            adj = [0 1 1 1; ...
                1 0 1 1; ...
                1 1 0 1; ...
                1 1 1 0];
        end

        function adj = graph_disconnected_components()
            % Create the graph:
            %     2         5
            %    / \       / \
            %   1---3     4---6 

            adj = [0 1 1 0 0 0; ...
                1 0 1 0 0 0; ...
                1 1 0 0 0 0; ...
                0 0 0 0 1 1; ...
                0 0 0 1 0 1; ...
                0 0 0 1 1 0];
        end

        function adj = graph_some_nodes_edgeless()
            % Create the graph:
            %          2     
            %         / \   
            %        4---3  
            %   
            %     Nodes 1, 5, 6 are edgeless

            adj = [0 0 0 0 0 0; ...
                   0 0 1 1 0 0; ...
                   0 1 0 1 0 0; ...
                   0 1 1 0 0 0; ...
                   0 0 0 0 0 0; ...
                   0 0 0 0 0 0];
        end

    end

    methods
        function verify_path_length(testCase, adjacency, startIdx, endIdx, expectedResult, debugTxt)
        
            % Execute the design 
            actualResult = shortest_path(adjacency, startIdx, endIdx);
        
            % Confirm the expected 
            msgTxt = sprintf('Failure context: %s', debugTxt);
            testCase.verifyEqual(actualResult, expectedResult, msgTxt);
        end
    end

end




