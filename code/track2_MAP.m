function [m x res] = track2_MAP(frames,bdryPts,ctrlPts,sizes,Sigma,Sigma_0,Sigma_1,dt,nofSteps,delta,nofIt,prior)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [m x res] = track2_MAP(frames,bdryPts,ctrlPts,sizes,Sigma,Sigma_0,Sigma_1,dt,nofSteps,delta,nofIt, prior)
%
% This function performs MAP estimation on each frame of the sequence. 
% It takes two contours, so it can be used for the heart images.
%
% Input:
%	frames - the sequence of images
%	bdryPts - the boundary points
%	ctrlPts - the control points
%	sizes - first element is the number of points in the outer circle
%	Sigma - the width of the noise kernel
%	Sigma_0 - the width of the deformation kernel
%	Sigma_1 - the width of the kernel of the ambient space
%	dt - time step of the ODE 
%	nofSteps - number of steps to generate the diffeomorphism
%       delta - time step in the optimization scheme
%	nofIt - number of iterations of gradient descent
%       prior - a binary variable which is 1 if we want to have a prior term and 0 for no prior (default = 1) 
%
% Output:
%
%	m - the final boundary points for each frame
%	x - the final control points for each frame 
%       res - a structure of the results
%
%
% Uses: MAP2.m, arcLength.m
%
% Note: in folder /project/track_all
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<12) 
  prior = 1;
end

% Setting input variables
res.input.frames = frames;
res.input.bdryPts = bdryPts;
res.input.ctrlPts = ctrlPts;
res.input.sizes = sizes;
res.input.Sigma = Sigma;
res.input.Sigma_0 = Sigma_0;
res.input.Sigma_1 = Sigma_1;
res.input.dt = dt;
res.input.nofSteps = nofSteps;
res.input.delta = delta;
res.input.nofIt = nofIt;
res.input.pror = prior;


% Some initial variables
nofBdryPts = size(bdryPts,1);
nofCtrlPts = size(ctrlPts,1);
nofFrames = size(frames(:,:,2:end),3);
E = zeros(nofIt,nofFrames);
mu = zeros(size(ctrlPts));
alpha = zeros(size(ctrlPts));
m = zeros(nofBdryPts,2,nofFrames);
x = zeros(nofCtrlPts,2,nofFrames);
m(:,:,1) = bdryPts;
x(:,:,1) = ctrlPts;


% calculating the likelihood paramaters
par_obs = parObsL2(frames(:,:,1),bdryPts,sizes,50,50,'off');
% par_obs(1:3,1)
% par_obs(1,1) = 0.20


figure(1)
colormap('gray')


for it = 1:nofFrames-1

  disp(it)  

  % reparameterizing by arc length
  if (1==0)
    step = 1/(sizes);
    [m1, nu1] = arcLength(m(1:sizes,:,it), step);  
    step = 1/(size(m,1)-sizes);
    [m2, nu2] = arcLength(m(sizes+1:end,:,it),step); 
    sizes(1,1) = size(m1,1);  
    m(:,:,it) = [m1;m2];
  end

  tic
  % finding the MAP
  [m(:,:,it+1) x(:,:,it+1) mu dummy dummy E(:,it)] = ...
    MAP2(frames(:,:,it+1),m(:,:,it),x(:,:,it),mu,alpha,sizes,...
    Sigma,Sigma_0,Sigma_1,dt,nofSteps,par_obs,delta,nofIt, prior);
  toc
  
  % displaying the result
  hold off
  imagesc(frames(:,:,it+1))
  hold on
  plot(m(:,1,it+1),m(:,2,it+1),'b.','Linewidth',2)
  plot(x(:,1,it+1),x(:,2,it+1),'ro','Linewidth',2)
  pause(0.01)
end


% storing the output
res.output.m = m;
res.output.x = x;
res.output.mu = mu;