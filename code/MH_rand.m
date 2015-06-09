function [m_new x_new w_new mu_new E_new w] = MH_rand(frame,bdryPts,ctrlPts,alpha,mu,Sigma,Sigma_0,Sigma_1,par_obs,dt,nofSteps,delta,nofIt,rnd,repar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [m_new x_new w_new mu_new E_new w] = MH(frame,bdryPts,ctrlPts,alpha,mu,Sigma,Sigma_0,Sigma_1,par_obs,dt,nofSteps,delta,nofIt,rnd)
%
% MH - this function implements the Metropolis-Hastings algorithm for sampling from the posterion distribution 
% The initial guess is alpha (in this case comes from the previous Importance sampling resampling procedure)
% 
% INPUT:
%	frame  - the image
%	bdryPts - the boundary points
%	ctrlPts - the control points
%	alpha - the initial guess for the distribution of alphas
%	mu - the mean for the prior on alpha
%	Sigma - 
%	Sigma_0 - 
%	Sigma_1 - 
%	par_obs
%	dt
%	nofSteps
%	nofIt
%
% OUTPUT:
%	m_new - the new positions of the boundary points
%	x_new - the new positions of the control points
%	w_new - the final state of the Markov chain
%	mu_new - the new transport of the mean to the new control points
%	w - the whole path of the Markov chain
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Some initialization
mu = zeros(size(mu));
w_new = zeros(size(mu));

nofCtrlPts = size(ctrlPts,1);
nofBdryPts = size(bdryPts,1);

% Constructing the covariance kernel matrix
C = K_matrix(ctrlPts,ctrlPts,2*Sigma_0^2 - 2*Sigma_1^2 + Sigma^2);
K = K_matrix(ctrlPts,ctrlPts,2*Sigma_0^2 - Sigma_1^2);

% Constructing the deformation kernel matrix
K_0 = K_matrix(ctrlPts,ctrlPts,Sigma_0^2);

% Constructing the covariance for the state
R = chol(C)*inv(K);

% Calculating the inverse of the covariance
Sigma_inv = K*inv(C)*K;

%mu = zeros(size(ctrlPts));
% Set the first state to be equal to the mean of the prior
w(:,:,1) = alpha;
m(:,:,1) = bdryPts;
x(:,:,1) = ctrlPts;



% Calculating the initial energy
  % calculating the prior term
  E_prior = (alpha(:) - mu(:))'*blkdiag(Sigma_inv,Sigma_inv)*(alpha(:)-mu(:))/2;
  % calulating the likelihood term
  E_obs = - logObsL_region(frame,bdryPts,par_obs,'off');
  E = E_prior + E_obs;
  E_min = E(1);

a = 0;

%  figure(1)
%  clf
%  imagesc(frame)
%  hold on
%  plot(bdryPts(:,1),bdryPts(:,2),'g','Linewidth',2)


mu_new = mu;
for it=1:nofIt
   
    % Generating a random state
    eta = delta*R'*rnd.states(:,:,it);
    w_temp = w(:,:,it) + eta;

    % Evolving along alpha
    [m_temp x_temp] = exponentialG(bdryPts,ctrlPts,w_temp,Sigma_0,dt,nofSteps);
    
    % calculating the prior term
    E_prior = (w_temp(:) - mu(:))'*blkdiag(Sigma_inv,Sigma_inv)*(w_temp(:)-mu(:))/2;
    % calulating the likelihood term
    E_obs = - logObsL_region(frame,m_temp,par_obs,'off');

    % calculating the total energy
    E_new = E_prior + E_obs;
    
    diff = E_new-E(it);
    u = rnd.r(it);
    rho = min(exp(-diff),1);
    if (u<rho)
      if(E_new<E_min)
      mode = m_temp;
      E_min = E_new;
      end
      a = a + 1; 
      E(it+1) = E_new;
      m(:,:,it+1) = m_temp;
      x(:,:,it+1) = x_temp;
      w(:,:,it+1) = w_temp;
      mu_new = w_new;
      %disp(E)
    else
      E(it+1) = E(it);
      m(:,:,it+1) = m(:,:,it);
      x(:,:,it+1) = x(:,:,it);
      w(:,:,it+1) = w(:,:,it);
    end

    m_new = m(:,:,it+1);
    x_new = x(:,:,it+1);
    w_new = w(:,:,it+1);
    %w_new = zeros(nofCtrlPts,2);
    
    

    % Calculating the acceptance rate
    rate = a/(it+1);

%      disp(rate)

    % Displaying the result
%      figure(1)
%      clf
 %     imagesc(frame)
    %  axis image
    %  hold on
    %  axis off
    %  plot(m(:,1,it), m(:,2,it), 'c', 'LineWidth', 2) ;
  %    plot(x(:,1,it), x(:,2,it), 'ro','Linewidth',2) ;
     % pause(0.01)
%      plot([mode(:,1);mode(1,1)], [mode(:,2); mode(1,2)], 'r', 'LineWidth', 2) ;

end 


if strcmp(repar,'on')
  % reparameterizing
  [m_new x_new] = reparam(m_new,x_new);
end


E_new = E(end);
mu_new = zeros(size(mu));
w_new = zeros(size(mu));
end