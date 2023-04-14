function [ik, ih, H_bst, K_bst] = get_hk_best(E00_all, ind_xi, H, K); 
%% Get the best H and K for a given index in E00_all 
% For run_hk_test_figs_merged.m 
E_tmp = E00_all{ind_xi }'; 
[E_max, imax] = max(E_tmp, [], 'all'); 
[ik, ih] = ind2sub(size(E_tmp), imax); 
H_bst = H(ih); 
K_bst = K(ik); 
end