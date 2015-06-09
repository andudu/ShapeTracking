


function [m_k_new x_k_new w_k_new alpha_k] = updateSample2(frame,m_k,x_k,mu_k,sizes,Sigma_small,Sigma_big,Sigma_0,Sigma_1,dt,nofSteps,par_obs,repar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% updateSample is a function which evolves the sample of particles along the exponential map, and calculates the
% the likelihood
% 
% INPUT: 
%	frame	 - the image
% 	m_k	 - the boundaries of the particles
%	x_k	 - the control points of the particles
%       mu_k	 - the mean of the state (it is currently set to zero from outside of the program)   	 
% 	Sigma_small - sigma for the overall covariance
%       Sigma_big - sigma for the block diagonal covariance
% 	Sigma_0  - sigma for projection	
%	Sigma_1  - sigma for shape space	
%	dt       - time step size
%	nofSteps - number of steps in ODE
%	par_obs  - obervationlikelihood parameters
% 	
% OUTPUT:
%	m_k_new  - the boundary pts after the reparameterization, before the evolution (TRICKY)
%       x_k_new  - the control pts after reparameterization before the evolution (TRICKY)
%	w_k_new	 - the weights of the particles
%       alpha_k - the alpha with which the particles get evolved
%
%	Warning: i am not really returning the deformed curves and control points!!!
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initializing some variables
nofBdryPts = size(m_k,1);
N = size(m_k,3);
nofCtrlPts = size(x_k,1);
m_k_new = zeros(nofBdryPts,2,N);
x_k_new = zeros(nofCtrlPts,2,N);
alpha_k= zeros(nofCtrlPts,2,N);
p = zeros(N,1);


disp('Updating the sample')


% generating the random numbers needed in the parloop
rnd = randn(nofCtrlPts,2,N);

% main parallel loop
parfor k = 1:N

  % calculating the covariance matrix
  K_0 = K_matrix(x_k(:,:,k),x_k(:,:,k),2*Sigma_0^2 - Sigma_1^2);
  % K = K_matrix(x_k(:,:,k),x_k(:,:,k),2*Sigma_0^2 - 2*Sigma_1^2 + Sigma^2 );
  % R = chol(K)*inv(K_0);
  
  if (~isempty(Sigma_big)) 
    % construct a non-homogeneous covariance matrix
    temp = x_k(:,:,k);
    K = Cov_block(temp(1:sizes(2),:),temp(sizes(2)+1:end,:),Sigma_small^2,Sigma_big^2);
  else 
    % construct a homogeneous covariance matrix
    K = K_matrix(x_k(:,:,k),x_k(:,:,k),2*Sigma_0^2 - 2*Sigma_1^2 + Sigma_small^2 );
  end
  R = chol(K)/K_0;
  



  % this needs to be fixed
  % projecting the velocity (think about the Sigmas)
  K01 = K_matrix(x_k(:,:,k),x_k(:,:,k),2*Sigma_0^2-Sigma_1^2);
  mu_k(:,:,k) = inv(K_0)*K01*mu_k(:,:,k);


  % generating random velocity
  alpha_k(:,:,k) = R'*rnd(:,:,k)+mu_k(:,:,k);
  
  % evolving along the group exponential
  [m_k_new(:,:,k) x_k_new(:,:,k)] = exponentialG(m_k(:,:,k),x_k(:,:,k),alpha_k(:,:,k),Sigma_0,dt,nofSteps);
  
  % evaluating the likelihood
  p(k,1) = logObsL_region2(frame,m_k_new(:,:,k),sizes,par_obs,'off');

  if strcmp(repar,'on')
    % reparameterizing
        temp1 = m_k_new(:,:,k);
        temp2 = x_k_new(:,:,k);
        [m1 x1] = reparam(temp1(1:sizes(1),:),temp2(1:sizes(2),:));
        [m2 x2] = reparam(temp1(sizes(1)+1:end,:),temp2(sizes(2)+1:end,:)); 
        m_k_new(:,:,k) = [m1;m2];
        x_k_new(:,:,k) = [x1;x2];
  end
end
% calculating the weights based on the observation likelihoods
w_k_new(:,1) = exp(p-max(p));

w_k_new = w_k_new/sum(w_k_new);
