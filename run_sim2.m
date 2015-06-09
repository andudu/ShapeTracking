%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script generates the simulations for the double contour example
% 

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
disp('1 - Double Gaussian' )
disp('2 - Double Salt & Pepper')
reply = input('');

if (reply~=0&&reply~=1&&reply~=2) 
  error('Your input should be an integer between 0 and 4.')
end


% creating a mat file to store the results
empty = [];
save 'results/simulations/paper_sim2' 'empty'


if (reply==0||reply==1)

% Experiment 1 - Gaussian noise 


% resetting the random stream
s = RandStream('mt19937ar','Seed',0)
RandStream.setGlobalStream(s)


% generating frames
[frames bdryPts ctrlPts alphas sim1] = sim_ellipse2(10,100,50,20,20,0.2,5,zeros(1,6),0.2);

size(bdryPts)
% tracking
[m_average x_average res1] = IS_Resample2(frames,bdryPts(:,:,1),ctrlPts(:,:,1),[100 10],20,[],20,20,0.2,5,1000);

% saving the first frame
imshow(frames(:,:,1),[min(min(min(frames))) max(max(max(frames)))],'border','tight')
print (gcf,'-dpng','results/simulations/double_gaussian0')
print (gcf,'-dps2','results/simulations/double_gaussian0')

% recording the results
track_video_MCMC2(res1.input.frames,res1.output.m_k,res1.output.x_k,100,'results/simulations/double_gaussian','results/simulations/double_gaussian');

save('results/simulations/paper_sim2','sim1','res1','-append')


end

if (reply==0||reply==2)

% Experiment 2 - Salt & pepper noise

% resetting the random stream
s = RandStream('mt19937ar','Seed',0)
RandStream.setGlobalStream(s)

% generating frames
[frames bdryPts ctrlPts alphas sim2] = sim_ellipse2(10,100,50,20,20,0.2,5,zeros(1,6),0.01);
frames = imnoise(frames,'salt & pepper',0.25);
sim2.output.frames = frames;
% tracking
[m_average x_average res2] = IS_Resample2(frames,bdryPts(:,:,1),ctrlPts(:,:,1),[100 10],20,[],20,20,0.2,5,1000);

% saving the first frame
imshow(frames(:,:,1),[min(min(min(frames))) max(max(max(frames)))],'border','tight')
print (gcf,'-dpng','results/simulations/double_sp0')
print (gcf,'-dps2','results/simulations/double_sp0')

% recording the results
track_video_MCMC2(res2.input.frames,res2.output.m_k,res2.output.x_k,100,'results/simulations/double_sp','results/simulations/double_sp');

save('results/simulations/paper_sim2','sim2','res2','-append')

end

rmpath(p1)
rmpath(p2)
matlabpool close
