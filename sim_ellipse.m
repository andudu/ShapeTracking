function [frames bdryPts ctrlPts alphas res] = sim_ellipse(nofCtrlPts,nofBdryPts,nofFrames,Sigma,Sigma_0,dt,nofSteps,par_aff,sig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [frames bdryPts ctrlPts alphas] = sim_ellipse(nofCtrlPts,nofBdryPts,nofFrames,Sigma,Sigma_0,dt,nofSteps,par_aff)
%
% This function generates a sequence of simulated shapes according to the following model:
%
% x_{t+1} = \phi_{t+1} ( A_{t+1} x_t )
% A_{t+1} = exp(noise)
%
% Both models are first order models.
% 
%
% 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% Storing the input variables
res.input.nofctrlPts = nofCtrlPts;
res.input.nofBdryPts = nofBdryPts;
res.input.nofFrames = nofFrames;
res.input.Sigma = Sigma;
res.input.Sigma_0 = Sigma_0;
res.input.dt = dt;
res.input.nofSteps = nofSteps;
res.input.par_aff = par_aff;

% storing the stream
defaultStream = RandStream.getDefaultStream;
savedState = defaultStream.State;
res.rand = savedState;

% Set parameters:

dim1 = 120; % the number of rows in the image
dim2 = 160; % the number of columnds in the image
R = 30; % the radius of the circle
center = [60 70];% the center of the circle

% Create some variables:
frames = zeros(dim1,dim2,nofFrames);
bdryPts = zeros(nofBdryPts,2,nofFrames);
nofCtrlPts = size((1:nofBdryPts/nofCtrlPts:nofBdryPts),2);
ctrlPts = zeros(nofCtrlPts,2,nofFrames);
alphas = zeros(nofCtrlPts,2,nofFrames);

% create an ellipse:
t = [0:2*pi/(nofBdryPts):2*pi]'; % might cause issues +-1 points
bdryPts(:,:,1) = [R*sin(t(1:end-1)) + center(1) (R-20)*cos(t(1:end-1))+center(2)]; % the last point will be repeated so we remove it
ctrlPts(:,:,1) = bdryPts(1:double(nofBdryPts)/(nofCtrlPts):end,:,1); % make sure control points are not repeated: important for the conditioning of the kernel


% nofCtrlPts = size(ctrlPts,1); % the number of control points might not be exact

BW = roipoly(frames(:,:,1), bdryPts(:,1,1),bdryPts(:,2,1));
frames(:,:,1) = 0.2*double(BW) + 0.3 + sig*randn(size(frames(:,:,1))); 
hold off
imagesc(frames(:,:,1))
axis off
axis tight
title(sprintf('frame %d',1))
hold on
colormap('gray')
plot(bdryPts(:,1,1),bdryPts(:,2,1),'.')
plot(ctrlPts(1:end,1,1),ctrlPts(1:end,2,1),'ro','Linewidth',2)
pause(0.1)

% generate motion


for i = 1:nofFrames-1

  
  [bdryPts(:,:,i+1),ctrlPts(:,:,i+1) ] = move_CP(bdryPts(:,:,i),ctrlPts(:,:,i),par_aff.*randn(1,6));
  
  % generate random vector field coefficients
  K_0 = K_matrix(ctrlPts(:,:,i+1),ctrlPts(:,:,i+1),Sigma_0^2);
  K_1 = K_matrix(ctrlPts(:,:,i+1),ctrlPts(:,:,i+1),Sigma^2);

  K_0_inv = inv(K_0);
  R = chol(K_1)*K_0_inv;
  disp(sprintf('The condition number of K_0 is %d',cond(K_0)))
  alphas(:,:,i) = R'*randn(nofCtrlPts,2); % the last one will be empty ( I can change that if I have a second order model ) 

  % Evolving the points
  [bdryPts(:,:,i+1) ctrlPts(:,:,i+1)] = exponentialG(bdryPts(:,:,i+1),ctrlPts(:,:,i+1),alphas(:,:,i),Sigma_0,dt,nofSteps);
  
  % Constructing the image
  frames(:,:,i+1) = roipoly(frames(:,:,i+1), bdryPts(:,1,i+1),bdryPts(:,2,i+1));
  frames(:,:,i+1) = 0.2*frames(:,:,i+1) +0.3 + sig*randn(size(frames(:,:,i+1))); 

  % Displaying the results
  hold off
  imagesc(frames(:,:,i+1))
  axis off
  axis tight
  title(sprintf('frame %d',i+1))
  hold on
  plot(bdryPts(:,1,i+1),bdryPts(:,2,i+1),'.-')
  plot(ctrlPts(:,1,i+1),ctrlPts(:,2,i+1),'ro','Linewidth',2)
  pause(0.1)

end


% Storing the output variables
res.output.frames = frames;
res.output.bdryPts = bdryPts;
res.output.ctrlPts = ctrlPts;
res.output.alphas = alphas;