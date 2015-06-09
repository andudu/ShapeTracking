function [m_k_new x_k_new alpha_k_new mu_k_new] = resample(m_k,x_k,alpha_k,mu_k,w_k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [m_k_new x_k_new w_k_new] = resample(m_k,x_k,w_k)
%
% a simple method for resampling the particles by sampling according to the 
% distribution of the weights
% unfortunately, this method does not allow parallelizing of the whole algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Resampling')
N = size(m_k,3);
m_k_new = zeros(size(m_k));
x_k_new = zeros(size(x_k));
alpha_k_new = zeros(size(alpha_k));
mu_k_new = zeros(size(mu_k));

c = cumsum(w_k);
u = zeros(N,1);
index = zeros(N,1);
parent = zeros(N,1);
u_tilde = rand/N;
i=1;
for j = 1:N
  %u_tilde = rand/N;
  u(j) = u_tilde+(j-1)/N;
  while (u(j)>c(i))
    i=i+1;
  end
  m_k_new(:,:,j) = m_k(:,:,i);%temp fix from i to i-1
  x_k_new(:,:,j) = x_k(:,:,i);
  alpha_k_new(:,:,j) = alpha_k(:,:,i);
  mu_k_new(:,:,j) = mu_k(:,:,i);
end
end