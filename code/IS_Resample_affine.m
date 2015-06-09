function m_k = IS_Resample_affine(frames,bdryPts,N,par_aff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [m_k x_k res] = IS_Resample_affine(frames, bdryPts,N, par_aff)
%
% This program implements Importance Sampling + Resampling with an affine motion prior.
% 
% INPUT:
%	frames - the image frames
%	bdryPts - the intitial boundary of the object
%	N - number of samples
%       par_aff - 6D variable containing the parameters of the affine transformation
%
% OUTPUT:
%	m_k - this is a (nofBdryPts,2,N,nofFrames) array which contains the sample of curves
%		representing the empirical estimate for the posterior at each time t
%
% USES: parObsL.m updateSample.m resample.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Some initial random variables
nofBdryPts = size(bdryPts,1);
nofFrames = size(frames(:,:,2:end),3);

% observation likelihood parameters
par_obs = parObsL(frames(:,:,1),bdryPts,100,100,'off');
par_obs(3,1) = 5;
par_obs(4,1) = 5;


% initializing the sample
m_k = zeros(nofBdryPts,2,N,nofFrames);
m_k(:,:,:,1) = reshape(repmat(bdryPts(:),1,N),nofBdryPts,2,N);


% initializing graphics
figure(1)
colormap('gray')
hold off


% Processing the frames
for it = 2:nofFrames
   
  disp(sprintf('Processing frame %d',it)) 
  
  tic
  % sample
  [m_k(:,:,:,it) w_k] = ...
	    updateSampleAffine(frames(:,:,it),m_k(:,:,:,it-1),par_obs,par_aff,'on');
  % resample
  [m_k(:,:,:,it)] = resampleAffine(m_k(:,:,:,it),w_k);
  toc

  % updating the center of mass
  par_obs(8) = mean(mean(squeeze(m_k(:,1,:,it))));
  par_obs(9) = mean(mean(squeeze(m_k(:,2,:,it))));

  % visualization
  hold off
  imagesc(frames(:,:,it))
  hold on
  plot(squeeze(m_k(:,1,:,it)),squeeze(m_k(:,2,:,it)))
  pause(0.01)

end


function m_k_new = resampleAffine(m_k,w_k)

disp('Resampling')
N = size(m_k,3);
m_k_new = zeros(size(m_k));

c = cumsum(w_k);
u = zeros(N,1);
u_tilde = rand/N;
i=1;
for j = 1:N
  u(j) = u_tilde+(j-1)/N;
  while (u(j)>c(i))
    i=i+1;
  end
  m_k_new(:,:,j) = m_k(:,:,i);
end


function [m_k_new w_k_new] = updateSampleAffine(frame,m_k,par_obs,par_aff,repar)

% Initializing some variables
nofBdryPts = size(m_k,1);
N = size(m_k,3);
m_k_new = zeros(nofBdryPts,2,N);
p = zeros(N,1);

disp('Updating the sample')

% main parallel loop
parfor k = 1:N
  % moving the points
  m_k_new(:,:,k) = move_affine(m_k(:,:,k),par_aff.*randn(1,6));
  % evaluating the likelihood
  p(k,1) = logObsL_region(frame,m_k_new(:,:,k),par_obs,'off');
end
   
% calculating the weights based on the observation likelihoods
p_max = max(p);
w_k = exp(p-p_max+100);
w_k_new = w_k/sum(w_k,1);


