function order = sortPatchOptions( sol, jMax, cell_ind )
% set the order to add to the queue
if size(jMax,1) > 1
    jMax = jMax';
end
acc_jMax = cumsum( [0 , jMax] );
chosen_idx = 1 + acc_jMax(cell_ind):acc_jMax(cell_ind + 1);
[~ , order] = sort( sol(chosen_idx) , 'descend');