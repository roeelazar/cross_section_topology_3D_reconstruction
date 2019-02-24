classdef BnBNode < handle
    % BnBNode Summary:
    %   holds all the parameters of a single node in the BnB search tree,
    %   e.i. holds all the parameters for a single relaxation we want to solve.
    %   Also, it holds the solution and its energy (LB), and the projected solution
    %   on the indicators space and its energy (UB).
    %
    % BnBNode Properties:
    %   ______Property___________________________________Description_____________________________________________________________________________type__________________
    %   | Root              | pointer to the root                                                                           | BnBroot 1 x 1                           |
    %   | Parent            | pointer the the node's parent                                                                 | BnBNode 1 x 1                           |
    %   | EQmatrix          | linear equalities constraints matrix: EQmatrix*x == EQvector                                  | double (|Cells| + 1) x num_of_vars      |
    %   | EQvector          | linear equalities constraints vector: EQmatrix*x == EQvector                                  | double (|Cells| + 1) x 1                |
    %   | INEQmatrix        | edge connectivity constraints INEQmatrix*x >= INEQvector  (num_of_vars x 1)                   | double edge constrtaints x num_of_vars  |
    %   | INEQvector        | edge connectivity constraints INEQmatrix*x >= INEQvector  (num_of_vars x 1)                   | double edge constrtaints x 1            |
    %   | Adj               | Adjacency matrix of each cell and pacth option                                                | cell |Cell| x 1                         |
    %   | Adj_size          | Size of the Adj matrix                                                                        | double 1 x 1                            |
    %   | Cost              | cost of each cell and patch option                                                            | double num_of_vars x 1                  |
    %   | Energy_bias       | energy bias resulting from fixed variables. must be added to the solver's output energy       | double 1 x 1                            |
    %   | RLX_UB            | energy of projected solution - **upper bound (Inf if the projected solution is infeasible)    | double 1 x 1                            |
    %   | X_UB              | projected solution on the indicators space (NaN if its infeasible)                            | double num_of_vars x 1                  |
    %   | RLX_LB            | energy of the relaxed solution - **lower bound (Inf if the solution is infeasible)            | double 1 x 1                            |
    %   | X_LB              | fractional relaxed solution (NaN if the its infeasible)                                       | double num_of_vars x 1                  |
    %   | Prune_branch      | boolean that indicate wether to prune the branch or not                                       | logical 1 x 1                           |
    %   | LP_iter_num       | counts the number of LP iteration                                                             | double 1 x 1                            |
    %   ---------------------------------------------------------------------------------------------------------------------------------------------------------------
    %  ** bounds refering to the optimization problem defined in BnBroot
    %
    % BnBNode Methods:
    %    BnBNode                 -  Constructor, seperates the cases of constructing a
    %                               root and branching a child.
    %    branchDown              -  Generating the children of a certain node with the 
    %                               identical properties of the parent beside 
    %                               Candidate property.
    %    dimensionReuction       -  Reducing the dimensions of the variables of the 
    %                               relaxation by fixing the patch of the last fixed cell,
    %                               and updatind all the optimization's parameters.
    %    findLeastCertainCellInd -  Among all unfixed cells return the cell with the 
    %                               maximal entropy. 
    %    solve_mincut            -  For the solution of the relaxation x, compute the 
    %                               partition of the vertices of the adjacency matrix 
    %                               induced by the min cut, from this partition return 
    %                               a linear inequality.
    %    projectSolution         -  In order to determine if the solution is an integer 
    %                               (as oppose to numerically integer), we project the
    %                               relaxed solution with round, and checking feasibiliy.
    %    expanded_x              -  The solution from the relaxation is dimensionally 
    %                               reduced, this method expand it to full length by
    %                               adding zeros and ones at the proper fixed cells slots.
    properties
        Root
        Parent
        Candidate
        EQmatrix
        EQvector
        INEQmatrix
        INEQvector
        Adj
        Adj_size
        Cost
        Energy_bias
        RLX_UB = Inf;
        X_UB = NaN;
        RLX_LB = Inf;
        X_LB = NaN;
        Prune_branch = false;
        LP_iter_num = 0;
    end
    
    methods
        function obj = BnBNode(Root, Parent)
            
            obj.Root = Root;
            if isempty(Parent)      % root node generation
                obj.Parent = [];
                obj.Candidate = nan(1, length(Root.jMax));
                obj.Parent.RLX_LB = -Inf;
                obj.EQmatrix = Root.EQmatrix;
                obj.EQvector = Root.EQvector;
                obj.INEQmatrix = [];
                obj.INEQvector = [];
                obj.Adj = Root.Adj;
                obj.Adj_size = Root.Adj_size;
                obj.Cost = Root.Cost;
                obj.Energy_bias = 0;
            else                    % non-root node generation
                obj.Parent = Parent;
                obj.Candidate = Parent.Candidate;
                obj.EQmatrix = Parent.EQmatrix;
                obj.EQvector = Parent.EQvector;
                obj.INEQmatrix = Parent.INEQmatrix;
                obj.INEQvector = Parent.INEQvector;
                obj.Adj = Parent.Adj;
                obj.Adj_size = Parent.Adj_size;
                obj.Cost = Parent.Cost;
                obj.Energy_bias = Parent.Energy_bias;
            end
            
        end
        
        function Children = branchDown(obj, cell_idx)
            % branching down from current node
            % note that the parameters are the same as the parent except
            % the candidate. the update occurs at the dimension reduction 
            Children = [];
            for patch_option = 1:obj.Root.jMax(cell_idx)
                child = BnBNode(obj.Root, obj);
                child.Candidate(cell_idx) = patch_option;
                Children = [Children, child];
            end
        end
                     
        function obj = dimensionReuction(obj)
            
            % the root doesn't need reduction
            if all(isnan(obj.Candidate))
                return
            end
            % get patch option of the last fixed cell
            undetermined_cells = isnan(obj.Candidate);
            Last_fixed_cell = setdiff(find(~isnan(obj.Candidate)), find(~isnan(obj.Parent.Candidate)));
            reduced_cell_idx = sum(undetermined_cells(1:Last_fixed_cell)) + 1;
            patch_option = obj.Candidate(Last_fixed_cell);
            
            % reduced_var_cell_idx indicate the location of the patch 
            % options of last fixed cell
            acc_jMax = cumsum([0 ; obj.Root.jMax(undetermined_cells)]);
            reduced_var_cell_idx = (acc_jMax(reduced_cell_idx)) + (1:obj.Root.jMax(Last_fixed_cell));            
            
            % --- reduce Cost and Energy_bias ---
            obj.Energy_bias = obj.Energy_bias + obj.Cost(reduced_var_cell_idx(patch_option));
            obj.Cost(reduced_var_cell_idx) = [];
            
            % --- reduce linear equality constraints ---
            
            % columns reduction
            
            % update vector
            obj.EQvector = obj.EQvector - obj.EQmatrix(:,reduced_var_cell_idx(patch_option));
            % remove matrix columns
            obj.EQmatrix(:,reduced_var_cell_idx) = [];
            
            % rows reduction
            
            zeroRows = ~any(obj.EQmatrix, 2);
            % necessary condition: if a row in the constraints matrix is
            % zeros, then the corresponding row in the updated vector must
            % be zero as well
            if any(obj.EQvector(zeroRows))
                obj.Prune_branch = true;
                return
            end
            % remove matrix rows
            obj.EQmatrix(zeroRows,:) = [];
            % remove vector rows
            obj.EQvector(zeroRows) = [];
            
            % --- reduce linear inequality constraints ---
            if ~isempty(obj.INEQvector)
                
                % columns reduction
                
                % ypdate vector
                obj.INEQvector = obj.INEQvector - obj.INEQmatrix(:,reduced_var_cell_idx(patch_option));
                % remove matrix columns
                obj.INEQmatrix(:,reduced_var_cell_idx) = [];
                
                % rows reduction
                
                % necessary condition: similar to the equality constraints
                zeroRows = ~any(obj.INEQmatrix, 2);
                if any(obj.INEQvector(zeroRows) > 0)
                    obj.Prune_branch = true;
                    return
                end
                % remove matrix rows
                obj.INEQmatrix(zeroRows,:) = [];
                % remove vector rows
                obj.INEQvector(zeroRows) = [];
            end
            
            % --- reduce Adj and update ConnCompSize ---
            obj = reduceAdj(obj, reduced_var_cell_idx, patch_option);
            
            
        end
        
        function idx = findLeastCertainCellInd(obj)
            % max entropy makes sense, since the varibles of each cell in
            % the relaxtion, are the probability simplex
            x = obj.X_LB;
            cell_matrix = obj.EQmatrix(1:length(obj.Root.jMax( isnan(obj.Candidate) )) , :);
            
            temp = -x.*log(x);
            temp(isnan(temp)) = 0;
            
            % max entropy
            [~,idx] = max( cell_matrix*temp );
            
            idx = find( cumsum( isnan(obj.Candidate))==idx );
        end
        
        function in_eq = solve_mincut(obj)         
            
            % construction of the adjacency of the solution x: Ax
            time2DuplicateX = cellfun( @(x) size(x,1) , obj.Adj);
            sol = obj.X_LB ; sol(sol < 0) = 0;  % entries (weights) of the adjacency must be non negative
            duplicatedX = cell2mat(  cellfun( @(x,y) x*ones(y,1) , num2cell( sol ) ,num2cell(time2DuplicateX) ,'UniformOutput' , false) );
            A_idx = cell2mat( obj.Adj );
            % sparse representation
            Ax = sparse( A_idx(:,1) , A_idx(:,2) , duplicatedX , obj.Adj_size , obj.Adj_size );
            % symmetrify the matrix
            Ax = Ax + Ax';
            
            % computing the partition induced by the min cut
            [ ~ , V2Remain] = mincutMex( full(Ax) );
            V2Remain = V2Remain';
            TwoSetsAdj = cellfun(@(x) V2Remain(x) ,obj.Adj ,'UniformOutput' , false);
            in_eq = cell2mat(cellfun(@xor2 , TwoSetsAdj ,'UniformOutput' , false) )';
            
        end
        
        function solve_sdp(obj)
            
            % generate Laplacian matrix
            n = obj.Adj_size;
            time2DuplicateX = cellfun( @(x) size(x,1) , obj.Adj);
            sol = obj.X_LB ; sol(sol < 0) = 0;
            duplicatedX = cell2mat(  cellfun( @(x,y) x*ones(y,1) , num2cell( sol ) ,num2cell(time2DuplicateX) ,'UniformOutput' , false) );
            A_idx = cell2mat( obj.Adj );

            Ax = sparse( A_idx(:,1) , A_idx(:,2) , duplicatedX , n , n );
            Ax = Ax + Ax';
            
            Laplacian = full( diag( Ax*ones(n, 1) ) - Ax );
            
            % check if the solution is already algebraicly connected
            fiedler_bound = 2*(1 - cos(pi/n));
            modified_Laplacian_eigvals = eig(Laplacian + fiedler/n);
            alg_connectecd = min(modified_Laplacian_eigvals)/fiedler_bound > 1;
            
            if alg_connectecd
                return
            end
            
            % solve sdp and update the lower bound
            yalmip('clear')
            num_of_vars = size( obj.EQmatrix, 1);
            
            % variables
            x_var = spdvar(num_of_vars, 1);
            
            % energy
            Energy = obj.Cost'*x_var;
            
            % constraints
            C = [ x_var(:)>=0 , obj.EQmatrix*x_var == obj.EQvector ,...
                obj.INEQmatrix*x_var >= obj.INEQvector, ...
                Laplacian >= fiedler_bound*(eye(n) - 1/n) ];
            
            % optimize - mosek sdp solver
            diagnostics = optimize(C, Energy, yalmip_ops);
            
            if ~any( ~([ 0 , 4 ] - diagnostics.problem))  % Problem with solver ( not numerical one )
                error(['Problem with solver: ',diagnostics.info])
            end
            
            obj.X_LB = double(x_var);
            obj.RLX_LB = double(Energy) + obj.Energy_bias;
        end
        
        function obj = projectSolution(obj, A, b, D, d)
            
            % intrinsic params
            feas_tol = 1e-5;
            
            % project to integer solution
            x_proj = round(obj.X_LB);
            
            % updating outputs
            obj.X_UB = x_proj;
            
            % check feasibility - equalities and inequalities constraints
            necessaryCond =  ( norm(A*x_proj - b,'inf') < feas_tol ) ;
            if ~isempty(D)
                necessaryCond = necessaryCond && ( sum(  0.5*( 1 - sign(D*x_proj - d + feas_tol)  ) ) <= 0.5 ) ;
            end
            
            % checking necessary condition before the expansive condition
            if ~necessaryCond
                return
            end
            
            % remain to check connectivity (expansive condition)
            A_idx = cell2mat( obj.Adj( logical( x_proj ) ) );
            
            % checking feasibity with number of connected components - no
            % need for expansive min cut
            G = graph( [A_idx(:,1); obj.Adj_size], [A_idx(:,2); obj.Adj_size]  );
            [~,binsize] = conncomp(G);
            
            if length(binsize) == 1 % condition: one connected component
                obj.RLX_UB = obj.Cost'*x_proj + obj.Energy_bias;
            end
            
            
        end
        
        function expanded_x = ExpandX(obj, reduced_x)

            jMax = obj.Root.jMax;
            candidate = obj.Candidate;
            
            if length(reduced_x)==sum(jMax) % solution is all ready at full length
                expanded_x = reduced_x ;
                return
            end
            
            expanded_x = zeros( sum(jMax) , 1);
            acc_jMax = cumsum([0 ; jMax]);
            determined_cells_idx = candidate + acc_jMax(1:end-1)';
            expanded_x(determined_cells_idx(~isnan(determined_cells_idx))) = 1;
            undetermined_cells_idx = logical( sum( obj.Root.EQmatrix(isnan(candidate),:) ) );
            expanded_x( undetermined_cells_idx ) = reduced_x ;
        end
        
    end
        
    
end

% Class related functions
    
    function out_obj = reduceAdj(obj, reduced_var_cell_idx, patch_option)
        % in order to reduce dimensions when we go down the tree we also
        % need to address the adjacency matrix. We take the fixed patch and
        % generating a new mapping of the current connected components of the
        % initial loops.
        
            out_obj = obj;
            connected_edges = out_obj.Adj{ reduced_var_cell_idx(patch_option) };
            out_obj.Adj( reduced_var_cell_idx ) = [];
            
            % Update only if a cell with connectivity is chosen
            if ~isempty(connected_edges)
            
                % Updating Adj - reindexation of vertices in Adj
                G = graph( [ connected_edges(:,1) ; max(max( cell2mat(out_obj.Adj))) ] ,[ connected_edges(:,2) ; max(max( cell2mat(out_obj.Adj))) ] );
                map = conncomp(G);
                out_obj.Adj = cellfun(@(x) map(x) , out_obj.Adj ,'UniformOutput' , false);
                
                % Updating Adj_size
                out_obj.Adj_size = max(map);
            
            end
            
            % checking feasibility of the reduced adjacency
            A_idx = cell2mat(out_obj.Adj);
            G = graph( [  A_idx(:,1) ; out_obj.Adj_size ] ,[  A_idx(:,2) ; out_obj.Adj_size ] );
            [~ , binsSizes] = conncomp(G);
            
            if length(binsSizes) > 1
                out_obj.Prune_branch = true;
                return
            end
            
        end

