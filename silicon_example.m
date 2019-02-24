
% INPUTS: 
%   ______Field_Name_________________________________Description_____________________________type____________
%   | V                 | list of points in R^3, sampled from the surface           | double |V| x 3        |
%   | E                 | list of edges that connect pair of vetrices from V        | double |E| x 2        |
%   | cellloop2edges    | defines the loops in each cell by list of edges from E    | cell   |Cells| x 1    |
%   | loopconnection    | dictates the connectivity of loops in each patch option   | cell   |Cells| x 1    |
%   |                   | refers to loops in cellloop2edges                         |                       |
%   | Cost              | the cost of each connectivity in loopconnection           | cell   |Cells| x 1    |
%   | genus             | target genus of the reconstructed surface                 | double 1 x 1          |
%   ---------------------------------------------------------------------------------------------------------

% OUTPUTS: 
%   ___relevant_feilds___________________description_____________________________________________type___________
%   | Sol               | the optimal, may be not optimal feasible solution             | double |Cells| x 1    |
%   |                   | if set FindOptimum = false, then you get the first feasible   |                       |
%   | solver.BnB_UB     | energy of Sol                                                 | double 1 x 1          |
%   | solver.Opt_time   | optimzation time                                              | double 1 x 1          |
%   --------------------------------------------------------------------------------------------------------------------------------------

%% load Data
load('silicium_data_example_planes57_genus217.mat')

%% preprocess
jMax = cellfun(@length, Cost);      % count the number of patch options per cell  
Cost = cell2mat(Cost);              % the cost stacked as a column vector
init_EC = size(V, 1) - size(E, 1);  % the Euler characteristics of the input undirected graph G(V,E) 
% contribution to Euler characteristics of each patch options
delta_EC = cell2mat(cellfun(@(z) cellfun( @(y) sum(cellfun(@(x) 2 - length(x), y )), z ), loopconnection, 'UniformOutput', false));

%% consruction of linear equalities constraints
acc_jMax = cumsum( [ 0 ; jMax]);
col_idx = 1:acc_jMax(end);
rows_idx = ones(1, length(col_idx));

for ii = 1:length(jMax)
    rows_idx((acc_jMax(ii)+1):acc_jMax(ii+1)) = ii ;
end

A = sparse( rows_idx, col_idx, ones(1,length(col_idx)), length(jMax), length(col_idx));

EQmatrix = [ A ; delta_EC'];
EQvector = [ ones(size(A,1) , 1 )   ; 2 - 2*genus - init_EC ];

%% construction of the adjacency matrix
% The data may not always be sampled from parallel planes, and the
% connected components may be more complex than simple loops structure.
% The connectivity should be between connected components.
[Adj, Adj_size] = get_Adj(E, cellloop2edges, loopconnection);

%% optimization
% setting parameters

% here you insert the solver you are using, for example I use mosek solver
% for other option checkout: https://yalmip.github.io/allsolvers/
% take into consideration that it should suppert mainly linear programming,
% and semidefinite programming (though without the SDP part, the alg is still
% robust)
YalmipOps = sdpsettings('solver','mosek','verbose',0,'cachesolvers',true);

Tmax = 10*60;           % set optimization time to 10 minutes
UB_thresh = 1e-4;       % numerical threshold for upper bound updating
FindOptimum = true;     % see above
Verbose = true;         % display info
solver = BnBroot(EQmatrix, EQvector, Adj, Adj_size, -Cost, jMax, UB_thresh, YalmipOps, Tmax, FindOptimum, Verbose);

BnBOptimization(solver)

if ~isfinite(solver.BnB_UB)
    disp('Infeasible solution')
    return
elseif FindOptimum && solver.Opt_time < Tmax
    disp('Optimal solution found')
else
    disp('Feasible solution found')
end
disp(['Optimization time: ',num2str(solver.Opt_time),' seconds'])
Sol = solver.X_opt2Sol;

%% SUBFUNCTIONS
function [Adj, Adj_size] = get_Adj(E, cellloop2edges, loopconnection)
% construction of the adjacency matrix 

G = graph(E(:, 1), E(:, 2));
V2conncomp = conncomp(G);
Adj_size = max(V2conncomp);

% reindexing: from local loops to connected compponent
loops2conncomp = cellfun(@(x) V2conncomp(E(cellfun(@(z) z(1) ,x ),1))' , cellloop2edges, 'UniformOutput', false);
connectivity_by_conncomp = cellfun(@(x,y) cellfun( @(z) cellfun(@(w) y(w), z, 'UniformOutput', false) , x, 'UniformOutput', false) , loopconnection, loops2conncomp, 'UniformOutput', false);
connectivity_by_conncomp = cellfun(@(x) cellfun(@(z) z(cellfun(@length, z) > 1), x, 'UniformOutput', false), connectivity_by_conncomp, 'UniformOutput', false); 

% generating the Adjacency matrix - storing indices
% since each pacth option represent combinatorial graph, storing indices
% is sufficient
Adj = {};
for ii = 1:length(connectivity_by_conncomp)
    Adj = [Adj ; connectivity_by_conncomp{ii} ];
end
Adj = cellfun(@(x) cell2mat(cellfun(@pair_Adj, x, 'UniformOutput', false)), Adj, 'UniformOutput', false);
end

function out = pair_Adj(x)
% arranging the Adj matrix by pair of indices
% if K components are connected, it will be represented as a path graph
% it will yeild tighter lower bound in the relaxation
if isempty(x) 
    out = x;
    return
end

out = [ x(1:(end-1)) , x(2:end) ];

end
