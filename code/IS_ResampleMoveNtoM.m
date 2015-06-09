function [m_k x_k res] = IS_ResampleMoveNtoM(frames, bdryPts,ctrlPts,Sigma,Sigma_0,Sigma_1,dt,nofSteps,N,M,delta,nofIt,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% [m_k x_k E_k res] = IS_ResampleMoveNtoM(frames, bdryPts,ctrlPts,Sigma,Sigma_0,Sigma_1,dt,nofSteps,N,M,nofIt,type)
%
% This program implements Importance Sampling + Resampling + MCMC move.
% 
% INPUT
%	frames 	- a sequence of images
%	bdryPts - bdryPts
%	ctrlPts - ctrlPts
%	Sigma 	- sigma for covariance
%	Sigma_0	- sigma for projected space
%	Sigma_1 - sigma for shape space
%	dt 	- time step in ODE
%	nofSteps - nofSteps in ODE
%	N 	- number of samples for IS
%       M 	- number of samples for MCMC
%	nofIt 	- number of iterations for MCMC moves
%	type   -  type of MCMC proposal distribution: RW - random walk; LD - Langevin Diffusion
%
% OUTPUT:
%	m_k - this is a (nofBdryPts,2,N,nofFrames) array which contains the sample of curves
%		representing the empirical estimate for the posterior at each time t
%	x_k - this is a (nofCtrlPts,2,N,nofFrames) array which contains the sample of ctrlPts
%      		at each time t
%       E_k - the energy for all particles
%       res - a structure containing all inputs, outputs, and the random stream
%       
%
% USES: parObsL.m updateSample.m resample.m MCMCmove_rand.m
%
% MCMCmove allows for the random stream to be stored
% MCMCmove.m can be called either with a MH-random-walk proposal or with a MH-Langevin proposal
% in this case - with MH_rand which accepts prespecified random streams
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%  global PAR
%  PAR.Sigma = Sigma;
%  PAR.Sigma_0 = Sigma_0;
%  PAR.Sigma_1 = Sigma_1;

% Setting input variables
res.input.frames = frames;
res.input.bdryPts = bdryPts;
res.input.ctrlPts = ctrlPts;
res.input.Sigma = Sigma;
res.input.Sigma_0 = Sigma_0;
res.input.Sigma_1 = Sigma_1;
res.input.dt = dt;
res.input.nofSteps = nofSteps;
res.input.N = N;
res.input.M = M;
res.input.delta = delta;
res.input.nofIt = nofIt;

% Storing the stream
defaultStream = RandStream.getGlobalStream;
savedState = defaultStream.State;
res.rand = savedState;

% Some initial random variables
nofBdryPts = size(bdryPts,1);
nofCtrlPts = size(ctrlPts,1);
nofFrames = size(frames(:,:,2:end),3);

par_obs = parObsL(frames(:,:,1),bdryPts,50,50,'off');

% initializing the sample
m_k = zeros(nofBdryPts,2,N,nofFrames);
x_k = zeros(nofCtrlPts,2,N,nofFrames);
E_k = zeros(N,nofFrames);
m_k(:,:,:,1) = reshape(repmat(bdryPts(:),1,N),nofBdryPts,2,N);
x_k(:,:,:,1) = reshape(repmat(ctrlPts(:),1,N),nofCtrlPts,2,N);


% initializing graphics
figure(1)
colormap('gray')
hold off

mu_k = zeros(nofCtrlPts,2,N);
mu_k = 0.1*randn(nofCtrlPts,2,N);

% Processing the frames
for it = 2:nofFrames
  
  disp(sprintf('Processing frame %d',it)) 
  mu_k = zeros(nofCtrlPts,2,N);
  
  tic
  % sample (the reparameterization in updateSample changes the location of bdryPts,ctrlPts)
  [m_k(:,:,:,it) x_k(:,:,:,it) w_k alpha_k] = ...
	    updateSample(frames(:,:,it),m_k(:,:,:,it-1),x_k(:,:,:,it-1),mu_k,Sigma,Sigma_0,Sigma_1,dt,nofSteps,par_obs,'off');
  

  % resample
  [m_k_dummy x_k_dummy alpha_k mu_k] = resampleNtoM(m_k(:,:,:,it-1),x_k(:,:,:,it-1),alpha_k,mu_k,w_k,M);

  % move
  [m_k_dummy x_k_dummy alpha_k mu_k E_k_temp] = ...
	    MCMCmove(frames(:,:,it),m_k_dummy,x_k_dummy,alpha_k,mu_k,Sigma,Sigma_0,Sigma_1,dt,nofSteps,par_obs,delta,nofIt,'on',type);


  % building blocks of the M particles to get to N particles needed for IS
  [m_k(:,:,:,it),x_k(:,:,:,it),alpha_k,mu_k, E_k(:,it)]= upsample(m_k_dummy,x_k_dummy,alpha_k,mu_k,N/M,E_k_temp);

  toc
  
  % visualization
  hold off
  imagesc(frames(:,:,it))
  hold on
  plot(squeeze(m_k(:,1,:,it)),squeeze(m_k(:,2,:,it)))
  plot(squeeze(x_k(:,1,:,it)),squeeze(x_k(:,2,:,it)),'o')
  pause(0.1)
end

res.output.m_k = m_k;
res.output.x_k = x_k;
res.output.E_k = E_k;
