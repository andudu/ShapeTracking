function [m_k_new w_k_new P_k_new]= update_affine_path(image,m_k,w_k,P_k,N,N_thr,edges,par_obs,par_aff)
nofBdryPts = size(m_k,1);
nofFrames = size(m_k,3);
m_k_new = zeros(nofBdryPts,2,N);
p = zeros(N,1);
figure(1)
hold off
imagesc(image)
axis off
hold on    
%disp('Generating new samples')
tic
parfor k=1:N
  %move 
    [m_k_new(:,:,k) A_k(:,:,k)] = move_affine_path(m_k(:,:,k),par_aff.*randn(1,6),'affine');

  %update weights
  %p(k,1) = p_obs_01_block(image,m_k_new(:,:,k),par_obs,'gaussian');
  
  p(k,1) = logObsL_region(image,m_k_new(:,:,k),par_obs,'block');

%    p(k) = p_obs_sc(image,m_k(:,:,k),p_max);
%  plot(m_k_new(:,1,k),m_k_new(:,2,k),'b','Linewidth',2)
  %plot(x_k_new(:,1,k),x_k_new(:,2,k),'ro','Linewidth',1)

end
toc

p_max = max(p);


%normalizing the weights
w_k_new(:,1) = w_k(:,1).*exp(p-p_max+100);
w_k_new = w_k_new/sum(w_k_new);

% resampling
[m_k_new w_k_new A_k_new P_k_new parent index] = resample(m_k_new,w_k_new,A_k,P_k);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Resampling function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m_k_new w_k_new A_k_new P_k_new parent index] = resample(m_k,w_k,A_k,P_k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [m_k_new x_k_new w_k_new] = resample(m_k,x_k,w_k)
%
% a simple method for resampling the particles by sampling according to the 
% distribution of the weights
% unfortunately, this method does not allow parallelizing of the algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(m_k,3);
m_k_new = zeros(size(m_k));

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
  index(j) = index(j)+1;
  parent(j) = i;
  A_k_new(:,:,j) = A_k(:,:,i);
  P_k_new(:,:,j) = A_k(:,:,i)*P_k(:,:,i);
end
w_k_new = ones(N,1)/N;
end
