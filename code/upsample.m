function [m_k_new,x_k_new,alpha_k_new,mu_k_new E_new]= upsample(m_k,x_k,alpha_k,mu_k,M,E)

% this function assumes we have a sample of size M and builds a sample of size kM
% by cooncatinating k blocks
% there should be way to do it by not expecting the new size to be a ,ultiple of the old size
% (i think I can just use resampling MtoN but the speed of the method will be slower)

% M - is the number of blocks

m_k_new = augment(m_k,M);
x_k_new = augment(x_k,M);
alpha_k_new = augment(alpha_k,M);
mu_k_new = augment(mu_k,M);
E_new = repmat(E,M,1);




function var_out = augment(var,M)
var1 = squeeze(var(:,1,:));
var2 = squeeze(var(:,2,:));
var1M = repmat(var1,1,M);
var2M = repmat(var2,1,M);
var_out(:,1,:) = var1M;
var_out(:,2,:) = var2M; 
