classdef BnBroot < handle
    % BnBroot Summary:
    %   this class holds all the parameters of the optimization problem we
    %   want to solve.
    %
    % BnBroot Properties:
    %    EQmatrix    - linear equalities matrix
    %    EQvector    - linear equalities vector
    %    Adj         - adjacency matrix (only indices)
    %    Adj_size    - size of the adjacency matrix
    %    Cost        - cost of each patch option
    %    jMax        - number of patch options per cell
    %    UB_thresh   - threshold of updating the upper bound (numerical purposes)
    %    YalmipOps   - yalmip options (solver, verbose, ...)
    %    Tmax        - max time for the optimization
    %    BnB_UB      - upper bound of the BnB method. Will be the energy of
    %                  the optimum at the end of the process. initialized to Inf
    %    X_opt       - the optimal solution
    %    FinfOptimum - boolean that indicates wether to go over all the tree, 
    %                  or to stop after the first feasible solution
    %    Verbose     - boolean, print the info
    %    Opt_time    - time of the optimization
    %    
    % BnBroot Methods:
    %    X_opt2Sol   - transforms the solution, from {0,1}^num_of_variable to
    %                  integer^|Cells|, i.e. transforms the solution of
    %                  indicators to an array of length |Cells| that contains 
    %                  integers, which represents the chosen patch option per cell
    
    properties
        EQmatrix
        EQvector
        Adj
        Adj_size
        Cost
        jMax
        UB_thresh
        YalmipOps
        Tmax
        BnB_UB = Inf;
        X_opt = NaN;
        FindOptimum;
        Verbose;
        Opt_time = 0;
    end
    
    methods
        % constructor
        function obj = BnBroot(EQmatrix, EQvector, Adj, Adj_size, Cost, jMax, UB_thresh, YalmipOps, Tmax, FindOptimum, Verbose)
            obj.EQmatrix = EQmatrix;
            obj.EQvector = EQvector;
            obj.Adj = Adj;
            obj.Adj_size = Adj_size;
            obj.Cost = Cost;
            obj.jMax = jMax;
            obj.UB_thresh = UB_thresh;
            obj.YalmipOps = YalmipOps;
            obj.Tmax = Tmax;
            obj.FindOptimum = FindOptimum;
            obj.Verbose = Verbose;
        end
        
        function Sol = X_opt2Sol(obj)
            Sol = zeros(1, length(obj.jMax));
            acc_jMax = cumsum([0; obj.jMax]);
            for ii = 1:length(obj.jMax) 
                Sol(ii) = find( obj.X_opt( (acc_jMax(ii) + 1):acc_jMax(ii + 1) ) );
            end
            
        end
    end
end

