function BnBOptimization(BnB_root)
% function solves the following optimization problem:
%---------------------------------------------------------------
%
% min_x             x^T * BnB_root.Cost
%
% such that:        x(:) in {0,1}
%                   BnB_root.EQmatrix * x == BnB_root.EQvector
%                   e(x) >= 1
%
%---------------------------------------------------------------
% where e(x) is the edge connectivity as defined in the paper
% input: BnB_root of type BnBroot

%% intrinsic params
counter = 0;
integer_solution_thresh = 1e-2;

%% Create root node of BnB tree
root = BnBNode(BnB_root, []) ;

% queue for the tree searching - BFS
optimization_que = CQueue(root);

optimization_start_time = tic;
while (optimization_que.size > 0) && (toc(optimization_start_time) <= BnB_root.Tmax)
    
    % pop the next node in the queue
    current_node = optimization_que.pop ;
        
    % necessary conditions to enter the relaxation
    if current_node.Parent.RLX_LB > ( BnB_root.BnB_UB + BnB_root.UB_thresh )
        continue
    end 
    
    % dimension reduction for current node (for relaxation)
    current_node.dimensionReuction;
    
    % necessary conditions resulting from the dimension reduction
    % dimensionReuction is expansive so this pruning occure afterward
    if current_node.Prune_branch
        continue
        
    % in the case we reach a leaf,
    % and ~Prune_branch, than the candidate is a feasible integer solution
    % its energy is Energy_bias
    % we update upper bound if necessary and continue the tree search
    elseif all(~isnan(current_node.Candidate)) && (current_node.Energy_bias < current_node.Root.BnB_UB - current_node.Root.UB_thresh)
        current_node.Root.BnB_UB = current_node.Energy_bias;
        current_node.Root.X_opt = current_node.ExpandX([]);
        continue
    end
    
    % updating counter (counting relaxations solved)
    counter = counter + 1;
    
    % solving relaxation
    solve_relaxation( current_node );
    
    % updating upper bound if nesessary
    integer_solution = (current_node.RLX_UB < Inf) && (norm(current_node.X_LB - round(current_node.X_LB),'inf') < integer_solution_thresh);

    if integer_solution && (current_node.RLX_UB < current_node.Root.BnB_UB - current_node.Root.UB_thresh)
        current_node.Root.BnB_UB = current_node.RLX_UB;
        current_node.Root.X_opt = current_node.ExpandX(current_node.X_UB);
    end
    
    % if integer feasible solution is found and we are only interested in
    % that, then return
    if (~current_node.Root.FindOptimum) && any(~isnan(current_node.Root.X_opt))
        return
    end 
  
    % printing the information of the iteration
    if BnB_root.Verbose
        fprintf('\n');
        fprintf('%4d.',counter);  
        fprintf('Fixed cells:    ');    fprintf('%4d',find(~isnan(current_node.Candidate)) ); fprintf('\n');
        fprintf('     Chosen patches: '); fprintf('%4d',current_node.Candidate(~isnan(current_node.Candidate)) ); fprintf('\n');
        fprintf('     LB =  %4.4f     UB = %4.4f\n',[current_node.RLX_LB , current_node.Root.BnB_UB]);
        fprintf('\n');
    end
    
    % --- Branching down ---    
    % condition for branching down:
    % 1.    The solution is fractional. If else, it was integer, and the
    %       optimal integer over the sub tree, then we update solution if
    %       needed and prune.
    % 2.    Did not reach a leaf. 
    % 3.    BnB condition - if the lower bound is less than the upper bound
    
    if ( ~integer_solution ) && any( isnan( current_node.Candidate ) ) && ( current_node.RLX_LB < ( BnB_root.BnB_UB + BnB_root.UB_thresh ) ) && (~current_node.Prune_branch)
        
        % find least certain cell (max entropy), add its "children" to the queue
        least_certain_cell_ind = current_node.findLeastCertainCellInd;
        
        Children = current_node.branchDown(least_certain_cell_ind);
        children_order = sortPatchOptions(current_node.ExpandX(current_node.X_LB), BnB_root.jMax, least_certain_cell_ind );
        for ii = 1:length(Children)
            optimization_que.push(Children(children_order(ii)));
        end
    end 
    
    % free memory
    current_node.Parent = [];
      
end

BnB_root.Opt_time = toc(optimization_start_time);

