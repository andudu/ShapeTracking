%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script generates the simulations for the single contour example 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% starting the matlab workers (no restriction on their number (up to 12))
m = matlabpool('size');
if (m==0)
  matlabpool open
end
 
p = pwd;
p1 = [p '/code'];
addpath (p1)
p2 = [p1 '/common'];
addpath (p2)



disp('Select an experiment:')
disp('0 - All')
disp('1 - Single Gaussian')
disp('2 - Single Gaussian with light occlusion')
disp('3 - Single Gaussian with dark occlusion' )
disp('4 - Single Salt & Pepper')
reply = input('');
if (reply~=0&&reply~=1&&reply~=2&&reply~=3&&reply~=4) 
  error('Your input should be an integer between 0 and 4.')
end


% creating a mat file to store the results
empty = [];
save 'results/simulations/paper_sim1' 'empty'


if (reply==0||reply==1)

% Experiment 1 - Gaussian noise 

disp('Gaussian noise')

s = RandStream('mt19937ar','Seed',1)
RandStream.setGlobalStream(s)


% generating frames
[frames bdryPts ctrlPts alphas sim1] = sim_ellipse(10,100,50,20,20,0.2,10,zeros(1,6),0.2);

% tracking
[m_average x_average res1] = IS_Resample(frames,bdryPts(:,:,1),ctrlPts(1:end-1,:,1),20,20,20,0.2,10,1000);


% saving the matlab variables 
save ('results/simulations/paper_sim1','sim1','res1','-append')

% saving the first frame

close(figure(1))
f = figure(1);
colormap('gray')
A = colormap('jet');
set(0,'DefaultAxesColorOrder',A)

hold off
imshow(frames(:,:,1),[min(min(min(frames))) max(max(max(frames)))],'border','tight')
print (gcf,'-dpng','results/simulations/single_gaussian0')
print (gcf,'-dps2','results/simulations/single_gaussian0')

% recording the results
track_video_MCMC(res1.input.frames,res1.output.m_k,res1.output.x_k,'results/simulations/single_gaussian','results/simulations/single_gaussian')


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 2 - Gaussian noise + white occlusion

if (reply==0||reply==2)
disp('Gaussian noise + light occlusion')


% resetting the random stream
s = RandStream('mt19937ar','Seed',1)
RandStream.setGlobalStream(s)


% generating frames
[frames bdryPts ctrlPts alphas sim2] = sim_ellipse(10,100,50,20,20,0.2,5,zeros(1,6),0.2);

x = 10;
y = 10;

for i = 1:50
  x = x + 1*(2*binornd(1,1)-1);
  y = y + 1*(2*binornd(1,1)-1);
  frames(y:y+15,x:x+15,i) = 0.5 + 0.2*randn(16,16);
  imagesc(frames(:,:,i))
  pause(0.01)
end

sim2.output.frames = frames;

% tracking
[m_average x_average res2] = IS_Resample(frames(:,:,1:end),bdryPts(:,:,1),ctrlPts(:,:,1),20,20,20,0.2,5,1000);


% saving the matlab variables 
save ('results/simulations/paper_sim1','sim2','res2','-append')

% saving the first frame
imshow(frames(:,:,30),[min(min(min(frames))) max(max(max(frames)))],'border','tight')
print (gcf,'-dpng','results/simulations/single_occluded_light0')
print (gcf,'-dps2','results/simulations/single_occluded_light0')

% recording the results
track_video_MCMC(res2.input.frames,res2.output.m_k,res2.output.x_k,'results/simulations/single_occluded_light','results/simulations/single_occluded_light')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 3 - Gaussian noise + dark occlusion

if (reply==0||reply==3)

disp('Gaussian noise + dark occlusion')

% resetting the random stream
s = RandStream('mt19937ar','Seed',1)
RandStream.setGlobalStream(s)


% generating frames
[frames bdryPts ctrlPts alphas sim3] = sim_ellipse(10,100,50,20,20,0.2,10,zeros(1,6),0.2);

x = 10;
y = 10;

for i = 1:50
  x = x + 1*(2*binornd(1,1)-1);
  y = y + 1*(2*binornd(1,1)-1);
  frames(y:y+15,x:x+15,i) = 0.2 + 0.2*randn(16,16);
  imagesc(frames(:,:,i))
  pause(0.01)
end

sim3.output.frames = frames;

% tracking
[m_average x_average res3] = IS_Resample(frames,bdryPts(:,:,1),ctrlPts(:,:,1),20,20,20,0.2,10,1000);

% saving the matlab variables 
save ('results/simulations/paper_sim1','sim3','res3','-append')


% saving the first frame
imshow(frames(:,:,30),[min(min(min(frames))) max(max(max(frames)))],'border','tight')
print (gcf,'-dpng','results/simulations/single_occluded_dark0')
print (gcf,'-dps2','results/simulations/single_occluded_dark0')

% recording the results
track_video_MCMC(res3.input.frames,res3.output.m_k,res3.output.x_k,'results/simulations/single_occluded_dark','results/simulations/single_occluded_dark')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 4 - Salt & pepper noise

if (reply==0||reply==4)
disp('Salt & peper noise')

% resetting the random stream
s = RandStream('mt19937ar','Seed',0)
RandStream.setGlobalStream(s)


% generating frames
[frames bdryPts ctrlPts alphas sim4] = sim_ellipse(10,100,50,20,20,0.2,5,zeros(1,6),0.01);
frames = imnoise(frames,'salt & pepper',0.25);
sim4.output.frames = frames;
% tracking
[m_average x_average res4] = IS_Resample(frames,bdryPts(:,:,1),ctrlPts(:,:,1),20,20,20,0.2,5,1000);

% saving the matlab variables 
save ('results/simulations/paper_sim1','sim4','res4','-append')


% saving the first frame
imshow(frames(:,:,1),[min(min(min(frames))) max(max(max(frames)))],'border','tight')
print (gcf,'-dpng','results/simulations/single_sp0')
print (gcf,'-dps2','results/simulations/single_sp')

% recording the results
track_video_MCMC(res4.input.frames,res4.output.m_k,res4.output.x_k,'results/simulations/single_sp','results/simulations/single_sp')

end

rmpath(p1)
rmpath(p2)
matlabpool close
