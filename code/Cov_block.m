 
function K = Cov_block(ctrlPts1,ctrlPts2,var_small, var_big)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% K = Cov_block(ctrlPts1,ctrlPts2,var_small, var_big)
%
% Uses:K_matrix.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K_small = K_matrix([ctrlPts1;ctrlPts2],[ctrlPts1;ctrlPts2], var_small);
K_big_block1 = K_matrix(ctrlPts1,ctrlPts1,var_big);
K_big_block2 = K_matrix(ctrlPts2,ctrlPts2,var_big);
K_big = blkdiag(K_big_block1,K_big_block2);

K = 0.5*K_small + 0.5*K_big;
