function [m_k x_k alpha_k mu_k E] =  MCMCmove(frame,m_k_old,x_k_old,alpha_k,mu_k,Sigma,Sigma_0,Sigma_1,dt,nofSteps,par_obs,delta,nofIt,repar,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% [m_k x_k alpha_k mu_k E] =  MCMCmove(frame,m_k_old,x_k_old,alpha_k,mu_k,Sigma,Sigma_0,Sigma_1,dt,nofSteps, par_obs,nofIt,repar,type)
%
% MCMCmove implements the Metropolis Hastings algorithm in parallel for each particle. 
% There are two possible implementations: one with random walk proposal and the other one with a Langevin diffusion proposal (selected through
% the type variable)
% 
% 
% Uses MH_rand.m MH_LD_rand.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



nofCtrlPts = size(x_k_old,1);
nofBdryPts = size(m_k_old,1);
N = size(m_k_old,3);
E = zeros(N,1);
disp('MCMCmove')

% mu_k = zeros(nofCtrlPts,2,N); % if this is uncommented: first order

% prespecifying the random variables I will use in Metropolis Hastings algorithm
rnd1 = randn(nofCtrlPts,2,nofIt,N);
rnd2 = rand(nofIt,N);
rnd = cell(N,1);
parfor k = 1:N
  rnd{k}.states = rnd1(:,:,:,k);
  rnd{k}.r = rnd2(:,k);

  if strcmp(type,'RW')
    % random walk proposal	
    [m_k(:,:,k) x_k(:,:,k) alpha_k(:,:,k) mu_k(:,:,k) E(k,1)] = ...
    MH_rand(frame,m_k_old(:,:,k),x_k_old(:,:,k),alpha_k(:,:,k),mu_k(:,:,k),Sigma,Sigma_0,Sigma_1,par_obs,dt,nofSteps,delta,nofIt,rnd{k},repar);
  elseif strcmp(type,'LD')
    % Langevin diffusion proposal
    [m_k(:,:,k) x_k(:,:,k) alpha_k(:,:,k) mu_k(:,:,k) E(k,1)] = ...
    MH_LD_rand(frame,m_k_old(:,:,k),x_k_old(:,:,k),alpha_k(:,:,k),mu_k(:,:,k),Sigma,Sigma_0,Sigma_1,par_obs,dt,nofSteps,delta,nofIt,rnd{k},repar);
  else 
    disp('Wrong type')
  end
end