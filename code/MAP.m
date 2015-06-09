function [m_new x_new w_new] = MAP(frame,bdryPts,ctrlPts,mu,Sigma,Sigma_0,Sigma_1,dt,nofSteps,par_obs,delta,nofIt,prior)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [m_new x_new w_new] = MAP(frame,bdryPts,ctrlPts,mu,Sigma,Sigma_0,Sigma_1,dt,nofSteps,par_obs,delta,nofIt,prior)
%
% This program finds the MAP estimate for a curve in an image, based on region based obs likelihood plus
% plus prior on the vector field.
% 
% Uses: arcLength.m, gradPhi.m, K_matrix.m
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note, I corrected the Sigma to Sigma_1 in the definition of K (11/05/2013)

% Some initial random variables
nofBdryPts = size(bdryPts,1);
nofCtrlPts = size(ctrlPts,1);

% Constructing the covariance kernel matrix
C = K_matrix(ctrlPts,ctrlPts,2*Sigma_0^2 - 2*Sigma_1^2 + Sigma^2);
K = K_matrix(ctrlPts,ctrlPts,2*Sigma_0^2 - Sigma_1^2);

% Constructing the deformation kernel matrix
K_0 = K_matrix(ctrlPts,ctrlPts,Sigma_0^2);

% Calculating the inverse of the covariance
Sigma_inv = K*(C\K);
Sigma_inv2 = blkdiag(Sigma_inv,Sigma_inv);


% initialize bdryPts and ctrlPts
m = zeros(nofBdryPts,2,nofIt);
x = zeros(nofCtrlPts,2,nofIt);
m(:,:,1) = bdryPts;
x(:,:,1) = ctrlPts;

% initial guess
alphas = zeros(nofCtrlPts,2);
mu = zeros(nofCtrlPts,2);

for it = 1:nofIt
  
  % calculating 
  [G m(:,:,it) x(:,:,it)] = gradPhi(bdryPts(:,:,1),ctrlPts(:,:,1),Sigma_0,dt,nofSteps,alphas);

  % calculate f at the new position
  f = - computeF(frame,m(:,:,it),par_obs);

  % parameterizing by arc length, and calculating the normal
  step = 1/size(m(:,:,it),1);
  [m(:,:,it), nu] = arcLength(m(:,:,it), step) ;

  % calculating the gradient 
  g = grad(m(:,:,it),nu,f);

  alphas = alphas - delta*reshape((G'*g(:) + prior*Sigma_inv2*(alphas(:) - mu(:))),nofCtrlPts,2);




 %   displaying the result
%    hold off
%    imagesc(frame)
%    hold on
%    plot(m(:,1,it),m(:,2,it),'b','Linewidth',2)
%    plot(x(:,1,it),x(:,2,it),'ro','Linewidth',2)
%    pause(0.01)
end
m_new = m(:,:,it);
x_new = x(:,:,it);
w_new = alphas;
end

function g = grad(m,nu,f)
 
    F = interp2(f, m(:,1), m(:,2), 'linear', max(max(f)));
    g(:,1) = (F.*nu(:,1)); 
    g(:,2) = (F.*nu(:,2));  
end

function f = computeF(frame,bdryPts,par_obs)
  
  % setting likelihood parameters
  cin = par_obs(1,1); cout = par_obs(2,1);
  varin = par_obs(3,1); varout = par_obs(4,1);
  
  
  %  cin = mean(V0(L==1));
  %  cout = mean(V0(L==0));
  

  % Calculating the energy inside and outside

  E_1 = (frame - cin).^2/(2*varin);
  E_2 = (frame - cout).^2/(2*varout);

  f = E_1 - E_2;
end

