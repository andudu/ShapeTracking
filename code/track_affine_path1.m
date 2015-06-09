 
function  [m_average w_k P_k] = track_affine_path1(frames,curves,N,N_thr,par_aff,movieName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%[m_average w_k P_k] = track_affine_path(frames,curves,N,N_thr,par_aff,movieName)
%
% INPUT:
%     frames - contains a sequence of curves
%     curves - the boundary points
%     N -  denotes the number of particles
%     N_thr - denotes the threshold number of particles
%     par_aff - the parameters of the affine transformation
%     movieName - the name of the movie
%
% OUTPUT: 
%         m_average - the mean curves
%         w_k - the eights if the affine transformations for each time
%         P_k - the affine transformation at each time
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Some initial variables
nofBdryPts = size(curves,1);
nofFrames = size(frames(:,:,2:end),3);
m_average = zeros(nofBdryPts,2,nofFrames);
w_k = zeros(N,nofFrames);



if (~isempty(movieName))
  mov = avifile([movieName,'.avi'], 'COMPRESSION', 'None', 'fps', 0.5, 'VideoName', movieName);
end

%Calculating the parameters of the observation density

par_obs = parObsL(frames(:,:,1),curves,100,100,'off');
par_obs(3,1) = 5;
par_obs(4,1) = 5;

% Find the edges in the video frames
edges = zeros(size(frames));
for i=1:nofFrames+1
  edges(:,:,i) = edge(frames(:,:,i),'canny',0.5);
end

%Displaying the first frame
figure(1)
colormap('gray')
imagesc(frames(:,:,2))
axis off
hold on
%Processing the first frame
disp(sprintf('Processing frame 1')) 
w_k(:,1) = ones(N,1)/N;
m = curves;

A_k = zeros(3,3,N);
P_k = zeros(3,3,N,nofFrames);
parfor k=1:N
  m_k(:,:,k) = m;
  [m_k(:,:,k) A_k(:,:,k)]= move_affine_path(m,par_aff.*randn(1,6),'affine');
%    plot(m_k(:,1,k),m_k(:,2,k),'b-','Linewidth',2)

  %p(k,1)=p_obs_01_block(frames(:,:,2),m_k(:,:,k),par_obs,'gaussian');
  p(k,1) = logObsL_region(frames(:,:,2),m_k(:,:,k),par_obs,'off');
  P_k(:,:,k,1) = eye(3);

end
	
p_max = max(p);
w_k(:,1) = w_k(:,1).*exp(p-p_max+100);
w_k(:,1) = w_k(:,1)/sum(w_k(:,1),1);


N_eff = 1/sum(w_k(:,1).^2);
disp(sprintf('N_eff= %d',(N_eff)))
%  N_thr = N*.75;%here I fix the threshold
%resampling if the effective sample size is too small

if (N_eff<N_thr)'
    disp('Many of the sample points are degenerate.')
    disp('Resampling...')
    [m_k w_k A_k P_k(:,:,:,1) parent index] = resample(m_k,w_k,A_k,P_k(:,:,:,1));
end


m_average(:,1,1) = squeeze(m_k(:,1,:))*w_k(:,1);
m_average(:,2,1) = squeeze(m_k(:,2,:))*w_k(:,1);

if (~isempty(movieName))
    disp('Saving the frame')
    frm = getframe(gcf);
    mov = addframe(mov,frm);
end

hold off
imagesc(frames(:,:,2))
axis off
hold on
plot(m_average(:,1),m_average(:,2),'r','Linewidth',2)

if (~isempty(movieName))
    disp('Saving the frame')
    frm = getframe(gcf);
    mov = addframe(mov,frm) ;
end

for i = 2:nofFrames
  tic
  disp(sprintf('Processing frame %d',i)) 

  %determining the new box parameters
  par_obs(8,1) = sum(m_average(:,1,i-1))/nofBdryPts;
  par_obs(9,1) = sum(m_average(:,2,i-1))/nofBdryPts;

  [m_k w_k(:,i) P_k(:,:,:,i)] = update_affine_path(frames(:,:,i+1),m_k,w_k(:,i-1),P_k(:,:,:,i-1),N,N_thr,edges(:,:,i+1),par_obs,par_aff);
  %[w I] = max(w_k(:,i));max is too sensitive to a particular curve
  %calculate the average
  disp(norm(m_average(:)))
  m_average(:,1,i) = squeeze(m_k(:,1,:))*w_k(:,i);
  m_average(:,2,i) = squeeze(m_k(:,2,:))*w_k(:,i);
  toc
   plot(m_average(:,1,i),m_average(:,2,i),'r','Linewidth',2)
  if (~isempty(movieName))
    disp('Saving the frame')
    frm = getframe(gcf);
    mov = addframe(mov,frm) ;
  end

    figure(1);colormap('gray')
    hold off
    imagesc(frames(:,:,i+1))
    axis off
    hold on
    plot(m_average(:,1,i),m_average(:,2,i),'r','Linewidth',2)
%      pause(0.5)

  if (~isempty(movieName))
    disp('Saving the frame')
    frm = getframe(gcf);
    mov = addframe(mov,frm) ;
  end
  
end
%  disp(w_k(:,end))
if (~isempty(movieName))
  mov = close(mov);
end
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
  A_k_new(:,:,j) = A_k(:,:,i);
  P_k_new(:,:,j) = A_k(:,:,i)*P_k(:,:,i);
  index(j) = index(j)+1;
  parent(j) = i;
end
w_k_new = ones(N,1)/N;
end
