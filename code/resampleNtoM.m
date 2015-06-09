function [m_k_new x_k_new alpha_k_new mu_k_new] = resampleNtoM(m_k,x_k,alpha_k,mu_k,w_k,M)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [m_k_new x_k_new alpha_k_new mu_k_new] = resampleNtoM(m_k,x_k,alpha_k,mu_k,w_k,M)
%
% This function performs resampling by changing the size of the new sample set from N to M
% 
% INPUT: the new sample size is M
% 	
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Resampling')

N = size(m_k,3);
nofBdryPts = size(m_k,1);
nofCtrlPts = size(x_k,1);
m_k_new = zeros(nofBdryPts,2,M);
x_k_new = zeros(nofCtrlPts,2,M);
alpha_k_new = zeros(nofCtrlPts,2,M);
mu_k_new = zeros(nofCtrlPts,2,M);
c = cumsum(w_k);
u = zeros(N,1);
index = zeros(N,1);
parent = zeros(N,1);
u_tilde = rand/N;
i=1;
for j = 1:M
  u(j) = u_tilde+(j-1)/M;
  while (u(j)>c(i))
    i=i+1;
  end
  m_k_new(:,:,j) = m_k(:,:,i); % temp fix from i to i-1
  x_k_new(:,:,j) = x_k(:,:,i);
  alpha_k_new(:,:,j) = alpha_k(:,:,i);
  mu_k_new(:,:,j) = mu_k(:,:,i);
end
end