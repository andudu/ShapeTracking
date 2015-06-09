


function [m_new x_new alpha_new mu_new alpha0_new E] = ...
MAP2(frame,bdryPts,ctrlPts,mu,alpha0,sizes,Sigma,Sigma_0,Sigma_1,dt,nofSteps,par_obs,delta,nofIt, prior)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 
% [m_new x_new alpha mu_new alpha0_new] = MAP2(frame,bdryPts,ctrlPts,mu,alpha0,sizes,Sigma,Sigma_0,Sigma_1,dt,nofSteps,par_obs,delta,nofIt,prior,repar)
%
% This program finds the MAP estimate for an object in an image, based on region based observation likelihood plus prior on the
% vector field. The object is representated by two closed contours and three regions: outside, inside and middle.
%
% Input:
%	frame - the image
%	bdryPts - the boundary points
%	ctrlPts - the control points
% 	mu - the mean of the prior
%       alpha0 - the starting point 
%	sizes - first element is the number of points in the outer circle
%	Sigma - the parameter of the noise kernel
%	Sigma_0 - the parameter of the deformation kernel
%	Sigma_1 - the parameter of the kernel of the ambient space
%	dt - time step
%	nofSteps - number of steps to generate the diffeomorphism
% 	par_obs - stores the mean intensity value for each region and their variances, ordered outside,middle, inside
%	nofIt - number of iterations of gradient descent
%       prior - a binary variable which is 1 if we want to have a prior term and 0 for no prior 
%
% Output:
%
%	m_new - new bdryPts
%	x_new - new ctrlPts
%       alpha_new - final solution
%       mu_new - (zeros!!!!)
%       alpha0_new - the new velocity (zeros)!!!!!
%       E - the final energy
%
% Uses:
%     gradPhi.m, logObs2_region2.m, K_matrix.m, exponentialG.m
%
% Note: in folder /project/track_all
%	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Some initial random variables
nofBdryPts = size(bdryPts,1);
nofCtrlPts = size(ctrlPts,1);
E = zeros(nofIt,1);

%par_obs = parObsL(frame,bdryPts,50,50,'off');

% Constructing the covariance kernel matrix
C = K_matrix(ctrlPts,ctrlPts,2*Sigma_0^2 - 2*Sigma_1^2 + Sigma^2);
K = K_matrix(ctrlPts,ctrlPts,2*Sigma_0^2 - Sigma_1^2);

% Constructing the deformation kernel matrix
K_0 = K_matrix(ctrlPts,ctrlPts,Sigma_0^2);
K_0_inv = inv(K_0);

% Calculating the inverse of the covariance
Sigma_inv = K*inv(C)*K;
m = zeros(nofBdryPts,2,nofIt);
x = zeros(nofCtrlPts,2,nofIt);

m(:,:,1) = bdryPts;
x(:,:,1) = ctrlPts;

%  figure(1)
%  colormap('gray')


% initial guess
alpha = alpha0;
% setting the mean
mu = zeros(nofCtrlPts,2);

for it = 1:nofIt

  % calculating 
  [G m(:,:,it) x(:,:,it)] = gradPhi(bdryPts(:,:,1),ctrlPts(:,:,1),Sigma_0,dt,nofSteps,alpha);

  % Calculate the energy
  L = logObsL_region2(frame,m(:,:,it),sizes,par_obs,'off');
  E(it) = - L + prior*alpha(:)'*blkdiag(Sigma_inv,Sigma_inv)*alpha(:)/2; 
  % disp(E(it))
  
  % calculate f at the new position
  [f_out f_in] = computeF1(frame,m(:,:,it),par_obs);

  % parameterizing by arc length, and calculating the normal
  % also reparameterizing the control points (this can cause trouble if a some second order model is used!) 
  step = 1/sizes;
  [m1, nu1] = arcLength(m(1:sizes,:,it), step);
  step = 1/(size(m,1)-sizes);
  [m2, nu2] = arcLength(m(sizes+1:end,:,it),step); 
  sizes(1,1) = size(m1,1);
  m(:,:,it) = [m1;m2];
  temp = m(1:(floor((nofBdryPts-1)/(nofCtrlPts-1))):end,:,it);
  if (size(temp,1)~=nofCtrlPts)
      x(:,:,it) = temp(1:nofCtrlPts,:,it);
  end
  x(:,:,it) = temp;

  % calculating the gradient 
  g = grad1(m1,m2,nu1,nu2,f_out,f_in);
  gg = reshape((G'*g(:) + prior*blkdiag(Sigma_inv,Sigma_inv)*(alpha(:) - mu(:))),nofCtrlPts,2);
 

% line search option, set ls=1 if needed
ls = 0;
if (ls==1)
    % linesearch
    E_new = E(it) + 1;
    d = delta;
    alpha_new = alpha;
    delta_min = 0.001;
    if(delta<delta_min)
      error('delta cannot be smaller than delta_min')
    end
    while (E_new>E(it)&&d>=delta_min)
      % update step
      alpha_new = alpha - d*gg;
      m_temp = exponentialG(bdryPts(:,:,1),ctrlPts(:,:,1),alpha_new,Sigma_0,dt,nofSteps);
      L = logObsL_region2(frame,m_temp,sizes,par_obs,'off');
      E_new = -L + prior*alpha_new(:)'*blkdiag(Sigma_inv,Sigma_inv)*alpha_new(:)/2;
      d = d/2;
    end
    alpha = alpha_new;
else 
    alpha = alpha - delta*gg;
end

end

%    % Displaying the result
%      hold off
%      imagesc(f_out)
%      hold on
%      plot(m(:,1,it),m(:,2,it),'b','Linewidth',2)
%      plot(x(:,1,it),x(:,2,it),'ro','Linewidth',2)
%      pause()

m_new = m(:,:,it);
x_new = x(:,:,it);
alpha_new = alpha;
mu_new = zeros(size(mu));
alpha0_new = zeros(size(alpha0));

end

% assumes there are only two regions: heart and background
function g = grad(m1,m2,nu1,nu2,f)
 
    F1 = interp2(f, m1(:,1), m1(:,2), 'linear', max(max(f)));
    F2 = interp2(f, m2(:,1), m2(:,2), 'linear', max(max(f)));

    g(:,1) = [F1.*nu1(:,1); - F2.*nu2(:,1)]; 
    g(:,2) = [F1.*nu1(:,2); - F2.*nu2(:,2)];

  
end
% assumes there are three regionsL: outside, ring, innner circle
function g = grad1(m1,m2,nu1,nu2,f_out,f_in)
 
    F1 = interp2(f_out, m1(:,1), m1(:,2), 'linear', max(max(f_out)));
    F2 = interp2(f_in, m2(:,1), m2(:,2), 'linear', max(max(f_in)));

    g(:,1) = [F1.*nu1(:,1);  F2.*nu2(:,1)]; 
    g(:,2) = [F1.*nu1(:,2);  F2.*nu2(:,2)];

  
end
% assumes there are only two regions: heart and background
function f = computeF(frame,bdryPts,par_obs)
% need to add block option and error checking option
%  
%  
% Label the image:  1 inside the contour; 0 outside the contour.
% L = roipoly(frame,bdryPts(:,1),bdryPts(:,2));
  
  
  % setting likelihood parameters
  cin = par_obs(2,1); cout = par_obs(1,1);
  varin = par_obs(4,1); varout = par_obs(3,1);
  
  
%  %  cin = mean(V0(L==1));
%  %  cout = mean(V0(L==0));
%  
%  % Calculating the energy inside and outside

  E_1 = (frame - cin).^2/(2*varin) + log(sqrt(2*pi*varin));
  E_2 = (frame - cout).^2/(2*varout) + log(sqrt(2*pi*varout));


f = E_1 - E_2;


end


function [f_out f_in] = computeF1(frame,bdryPts,par_obs)
  % function [f_out f_in] = computeF1(frame,bdryPts,par_obs)
  % cumputeF1 calculates grad f (f is the function inside the integral which needs to be minimized)
  % it assumes there are three regions: outside, ring, inner circle
  % (need to add block option and error checking option)
  
  
  % setting likelihood parameters
  mu1 = par_obs(1,1); mu2 = par_obs(2,1); mu3 = par_obs(3,1);
  var1 = par_obs(4,1); var2 = par_obs(5,1); var3 = par_obs(6,1);

  %  used only if we allow for the means of the region to change 
  %  cin = mean(V0(L==1));
  %  cout = mean(V0(L==0));

  
  % Calculating the energy inside and outside
  E_1 = (frame - mu1).^2/(2*var1) + log(sqrt(2*pi*var1));
  E_2 = (frame - mu2).^2/(2*var2) + log(sqrt(2*pi*var2));
  E_3 = (frame - mu3).^2/(2*var3) + log(sqrt(2*pi*var3));

  f_out = E_2 - E_1;
  f_in = E_3 - E_2;

end
