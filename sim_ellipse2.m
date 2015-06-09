function [frames bdryPts ctrlPts alphas res] = sim_ellipse2(nofCtrlPts1,nofBdryPts1,nofFrames,Sigma,Sigma_0,dt,nofSteps,par_aff,sig)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [frames bdryPts ctrlPts alphas] = sim_ellipse2(nofCtrlPts1,nofBdryPts1,nofFrames,Sigma,Sigma_0,dt,nofSteps,par_aff)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% Storing the input variables
res.input.nofctrlPts1 = nofCtrlPts1;
res.input.nofBdryPts1 = nofBdryPts1;
res.input.nofFrames = nofFrames;
res.input.Sigma = Sigma;
res.input.Sigma_0 = Sigma_0;
res.input.dt = dt;
res.input.nofSteps = nofSteps;
res.input.par_aff = par_aff;
res.input.sig = sig;

% storing the stream
defaultStream = RandStream.getDefaultStream;
savedState = defaultStream.State;
res.rand = savedState;

% Set parameters:

dim1 = 120; % the number of rows in the image
dim2 = 160; % the number of columnds in the image
R = 30; % the radius of the circle
center = [50 70];% the center of the circle

% inner contour variables
nofBdryPts2 = 75; % these need to be defined by the user
nofCtrlPts2 = 6; 
bdryPts2 = zeros(nofBdryPts2,2,nofFrames);
nofCtrlPts2 = size((1:nofBdryPts2/nofCtrlPts2:nofBdryPts2),2);
ctrlPts2 = zeros(nofCtrlPts2,2,nofFrames);

% Initialize variables:
frames = zeros(dim1,dim2,nofFrames);
bdryPts1 = zeros(nofBdryPts1,2,nofFrames);
nofCtrlPts1 = size((1:nofBdryPts1/nofCtrlPts1:nofBdryPts1),2);
ctrlPts1 = zeros(nofCtrlPts1,2,nofFrames);
alphas = zeros(nofCtrlPts1+nofCtrlPts2,2,nofFrames);

nofBdryPts = nofBdryPts1 + nofBdryPts2;
nofCtrlPts = nofCtrlPts1 + nofCtrlPts2;
ctrlPts = zeros(nofCtrlPts,2,nofFrames);
bdryPts = zeros(nofBdryPts,2,nofFrames);

% create first ellipse:
t = [0:2*pi/(nofBdryPts1):2*pi]'; % might cause issues +-1 points
bdryPts1(:,:,1) = [R*sin(t(1:end-1)) + center(1) (R-10)*cos(t(1:end-1))+center(2)]; % the last point will be repeated so we remove it
ctrlPts1(:,:,1) = bdryPts1(1:double(nofBdryPts1)/(nofCtrlPts1):end-1,:,1); % make sure control points are not repeated: important for the conditioning of the kernel

% create second ellipse:
t = [0:2*pi/(nofBdryPts2):2*pi]';
bdryPts2 = [20*sin(t(1:end-1)) + center(1) (15)*cos(t(1:end-1))+center(2)];
ctrlPts2(:,:,1) = bdryPts2(1:double(nofBdryPts2)/nofCtrlPts2:end-1,:,1);



% nofCtrlPts = size(ctrlPts,1); % the number of control points might not be exact

% Constructing the image
frames(:,:,1) = color_image(frames(:,:,1),bdryPts1(:,:,1),bdryPts2(:,:,1));


figure(1)
hold off
imagesc(frames(:,:,1))
axis off
axis tight
title(sprintf('frame %d',1))
hold on
colormap('gray')
bdryPts(:,:,1) = [bdryPts1(:,:,1);bdryPts2(:,:,1)];
ctrlPts(:,:,1) = [ctrlPts1(:,:,1);ctrlPts2(:,:,1)];
plot(bdryPts1(:,1,1),bdryPts1(:,2,1),'.')
plot(bdryPts2(:,1,1),bdryPts2(:,2,1),'.')
plot(ctrlPts1(:,1,1),ctrlPts1(:,2,1),'ro','Linewidth',2)
plot(ctrlPts2(:,1,1),ctrlPts2(:,2,1),'ro','Linewidth',2)

% generate motion


for i = 1:nofFrames-1
  
  [bdryPts(:,:,i+1),ctrlPts(:,:,i+1)] = move_CP(bdryPts(:,:,i),ctrlPts(:,:,i),par_aff.*randn(1,6));
  
  % generate random vector field coefficients
  K_0 = K_matrix(ctrlPts(:,:,i+1),ctrlPts(:,:,i+1),Sigma_0^2);
  K_1 = K_matrix(ctrlPts(:,:,i+1),ctrlPts(:,:,i+1),Sigma^2);

  K_0_inv = inv(K_0);
  R = chol(K_1)*K_0_inv;
  disp(sprintf('The condition number of K_0 is %d',cond(K_0)))
  alphas(:,:,i) = R'*randn(nofCtrlPts,2); % the last one will be empty ( I can change that if I have a second order model ) 

  % Evolving the points
  [bdryPts(:,:,i+1) ctrlPts(:,:,i+1)] = exponentialG(bdryPts(:,:,i+1),ctrlPts(:,:,i+1),alphas(:,:,i),Sigma_0,dt,nofSteps);

  bdryPts1(:,:,i+1) = bdryPts(1:nofBdryPts-nofBdryPts2,:,i+1);
  bdryPts2(:,:,i+1) = bdryPts(nofBdryPts-nofBdryPts2+1:end,:,i+1);  
  ctrlPts1(:,:,i+1) = ctrlPts(1:nofCtrlPts-nofCtrlPts2,:,i+1);
  ctrlPts2(:,:,i+1) = ctrlPts(nofCtrlPts-nofCtrlPts2+1:end,:,i+1);
  
  % Constructing the image
  frames(:,:,i+1) = color_image(frames(:,:,i+1),bdryPts1(:,:,i+1),bdryPts2(:,:,i+1));
  %frames(:,:,i+1) = imnoise(frames(:,:,i+1),'salt & pepper',0.2);
  % Displaying the results
  hold off
  imagesc(frames(:,:,i+1))
  title(sprintf('frame %d',i))
  hold on
  plot(bdryPts(:,1,i+1),bdryPts(:,2,i+1),'.')
  plot(ctrlPts(:,1,i+1),ctrlPts(:,2,i+1),'ro','Linewidth',2)
  pause(0.01)
end


% Add noise to the frames
% frames = frames + 0.5*randn;
%frames = imnoise(frames,'salt&pepper',0.5);

% Storing the output variables
res.output.frames = frames;
res.output.bdryPts = bdryPts;
res.output.ctrlPts = ctrlPts;
res.output.alphas = alphas;


function im_out = color_image(im_in,bdryPts,bdryPts1)
[dummy BW] = roifill(im_in,bdryPts(:,1),bdryPts(:,2));
im_out = 0.5*BW+(1-BW)*0.2;
[dummy BW] = roifill(im_out, bdryPts1(:,1),bdryPts1(:,2));
im_out = (1-BW).*im_out + BW*0.7;
im_out = im_out + sig*randn(size(im_out));
end

end
