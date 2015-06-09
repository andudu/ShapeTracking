function [m_k_new x_k_new w_k_new alpha_k parent] = update_CP_sd(image,m_k,x_k,w_k,N,N_thr,Sigma,Sigma_0,dt,nofSteps,par_obs,edges,P,mov)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [m_k_new x_k_new w_k_new alpha_k parent] = update_CP(image,m_k,x_k,w_k,N,N_thr,Sigma,Sigma_0,dt,nofSteps,par_obs,edges,P,mov)
%
% This function first moves the shape using the supplied affine transformation in P, and then deforms the shape through a deffeomorphic transformaiton.
%
% INPUT: image - the current frame (used only for visualization)
%        m_k - the current boundary points
%        x_k - the current control points
%        w_k - the current weights
%        N - the number of particles
%        N_thr - the threshold number of particles
%        Sigma - the sigma of the noise kernel
%        Sigma_0 - the sigma of the deformation kernel
%        dt - the timestep
%        nofSteps - number of small deformations
%        par_obs - parameters of the observation likelihood
%        edges - the edges in the image (needed for the PPP observation likelihood)
%        P - the affine transformation at time t
%        mov - the name of the movie
%    
% OUTPUT:m_k_new - the new boundary points
%        x_k_new - the new control points
%        w_k_new - the new weights
%        alpha_k - the coefficients of the deformations
%        parent - the index of the parent particle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nofBdryPts = size(m_k,1);
nofFrames = size(m_k,3);
nofCtrlPts = size(x_k,1);
m_k_new = zeros(nofBdryPts,2,N);
x_k_new = zeros(nofCtrlPts,2,N);
p = zeros(N,1);
parent = [1:N]';

figure(1)
hold off
imagesc(image)
axis off
hold on  
  
%disp('Generating new samples')
alpha_k = zeros(nofCtrlPts,2,N);
parfor k=1:N
  %Constructing the matrix of covariances to generate the random alphas
  x_k_temp = x_k(:,:,k);
  z2 = x_k_temp(:,1).^2*ones(1,nofCtrlPts);
  dist =  z2 - 2* x_k_temp(:,1)*x_k_temp(:,1)' + z2' ;
  z2 = x_k_temp(:,2).^2 * ones(1,nofCtrlPts) ;
  dist =  dist + z2 - 2* x_k_temp(:,2)*x_k_temp(:,2)' + z2' ;
  K_0 = exp(-dist/(2*Sigma_0*Sigma_0))/(2*pi*Sigma_0*Sigma_0);
  K_1 = exp(-dist/(2*Sigma*Sigma))/(2*pi*Sigma*Sigma);
  K_0_inv = inv(K_0);
  R = chol(K_1)*K_0_inv;
%    disp(sprintf('The condition number of K_0 is %d',cond(K_0)))
  alpha_k = reshape(R'*randn(nofCtrlPts,2*N),nofCtrlPts,2,N);

  %move 
  [m_k_new(:,:,k) x_k(:,:,k)] = move_CP_sd(m_k(:,:,k),x_k(:,:,k),P);
  
  %deform
  x_k_new(:,:,k) = x_k(:,:,k);
  for step = 1:nofSteps   
      [m_k_new(:,:,k) x_k_new(:,:,k)]= deform_CP(m_k_new(:,:,k),x_k(:,:,k),x_k_new(:,:,k),alpha_k(:,:,k),Sigma,Sigma_0,dt);
  end
  x_k(:,:,k)  = x_k_new(:,:,k);



%    
 % imagesc(image)
%    title(strcat('Sigma\_cp=',sprintf('%d, ',Sigma_0),'Sigma\_ns=',sprintf('%d, ',Sigma),'dt=',sprintf('%d, ',dt),'nofSteps=',sprintf('%d',nofSteps)))
%      plot(m_k_new(:,1,k),m_k_new(:,2,k),'b-','Linewidth',2)
%      plot(x_k_new(:,1,k),x_k_new(:,2,k),'ro','Linewidth',2)
  %pause(0.01)
 
  p(k,1) = logObsL_region(image,m_k_new(:,:,k),par_obs,'block');

   
end

p_max = max(p);
w_k_new(:,1) = w_k(:,1).*exp(p-p_max+100);


% normalizing the weights
w_k_new = w_k_new/sum(w_k_new);
% resampling
[m_k_new x_k_new w_k_new parent index] = resample(m_k_new,x_k_new,w_k_new);


function [m_k_new x_k_new w_k_new parent index] = resample(m_k,x_k,w_k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [m_k_new x_k_new w_k_new] = resample(m_k,x_k,w_k)
%
% a simple method for resampling the particles by sampling according to the 
% distribution of the weights
% unfortunately, this method does not allow parallelizing of the whole algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = size(m_k,3);
m_k_new = zeros(size(m_k));
x_k_new = zeros(size(x_k));
%first construct the CDF
%  c = zeros(N,1);
%  c(1) = w_k(1);
%  for i=2:N
%    c(i) = c(i-1)+w_k(i);
%  end
c = cumsum(w_k);
u = zeros(N,1);
index = zeros(N,1);
parent = zeros(N,1);
u_tilde = rand/N;
i=1;
for j = 1:N
  u(j) = u_tilde+(j-1)/N;
  while (u(j)>c(i))
    i=i+1;
  end
  m_k_new(:,:,j) = m_k(:,:,i);%temp fix from i to i-1
  x_k_new(:,:,j) = x_k(:,:,i);
  index(j) = index(j)+1;
  parent(j) = i;
end
w_k_new = ones(N,1)/N;
end

end