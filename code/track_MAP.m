function [m x res] = track_MAP(frames,bdryPts,ctrlPts,Sigma,Sigma_0,Sigma_1,dt,nofSteps,delta,nofIt,prior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [m x res] = track_MAP(frames,bdryPts,ctrlPts,Sigma,Sigma_0,Sigma_1,dt,nofSteps,delta,nofIt,prior)
% 
% This program implements a sequential segmentation algorithm on a sequence of images. It aims to maximize
% the posterior at each step, assuming at the previous step a correct maximizer was found.
%
% INPUT:
%	frames    - the sequence of images
%	bdryPts   - the boundary points
%	ctrlPts   - the control points
%	sizes     - first element is the number of points in the outer circle
%	Sigma     - the width of the noise kernel
%	Sigma_0   - the width of the deformation kernel
%	Sigma_1   - the width of the kernel of the ambient space
%	dt        - time step of the ODE 
%	nofSteps  - number of steps to generate the diffeomorphism
%       delta     - time step in the optimization scheme
%	nofIt     - number of iterations of gradient descent
%       prior     - a binary variable: 1 if we want to have a prior term and 0 for no prior (default = 1) 
%
% OUTPUT:
%	m   - the final boundary points for each frame
%	x   - the final control points for each frame 
%	res - contains the inputs and the outputs of the algorithm
%
% USES:
%	MAP.m, par_ObsL.m, reparam.m
%
%
% Note: this program returns the bdryPts and the ctrlPts after reparameterization 
%      (thus not the final contour of the exponential map (but the initial contour for the next step))
%      in folder /project/track_all
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<11) 
  prior = 1;
end


% Setting input variables
res.input.frames = frames;
res.input.bdryPts = bdryPts;
res.input.ctrlPts = ctrlPts;
res.input.Sigma = Sigma;
res.input.Sigma_0 = Sigma_0;
res.input.Sigma_1 = Sigma_1;
res.input.dt = dt;
res.input.nofSteps = nofSteps;
res.input.delta = delta;
res.input.nofIt = nofIt;
res.input.nofIt = prior;

% Some initial variables
nofBdryPts = size(bdryPts,1);
nofCtrlPts = size(ctrlPts,1);
nofFrames = size(frames(:,:,2:end),3);
mu = zeros(nofCtrlPts,2,nofFrames);
m = zeros(nofBdryPts,2,nofFrames);
x = zeros(nofCtrlPts,2,nofFrames);
m(:,:,1) = bdryPts;
x(:,:,1) = ctrlPts;

% calculating the likelihood paramaters
par_obs = parObsL(frames(:,:,1),bdryPts,50,50,'off');

figure(1)
colormap('gray')

for it = 1:nofFrames

  disp(sprintf('Processing frame %d',it))  
  
  % reparameterizing the curve
  [m(:,:,it) x(:,:,it)] = reparam(m(:,:,it),x(:,:,it));
 
  %calculating the maximum aposteriori estimate at each step
  [m(:,:,it+1) x(:,:,it+1) mu(:,:,it)] = ...
    MAP(frames(:,:,it+1),m(:,:,it),x(:,:,it),mu,Sigma,Sigma_0,Sigma_1,dt,nofSteps,par_obs,delta,nofIt,prior);
  
  
  % displaying the result
  hold off
  imagesc(frames(:,:,it+1))
  hold on
  plot(m(:,1,it+1),m(:,2,it+1),'b','Linewidth',2)
  plot(x(:,1,it+1),x(:,2,it+1),'ro','Linewidth',2)
  pause(0.01)

  % storing the output
  res.output.m = m;
  res.output.x = x;
  res.output.mu = mu;

end