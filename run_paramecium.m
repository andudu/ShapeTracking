% starting the matlab workers (no restriction on their number (up to 12))
m = matlabpool('size');
if (m==0)
  matlabpool open
end
% setting the random stream
s = RandStream('mt19937ar','Seed',0)
RandStream.setGlobalStream(s)

% loading the data
cd data
disp('Loading the data')
load('param_noisy.mat')

cd ..
% run the experiments 
p = pwd;
p1 = [p '/code'];
addpath (p1)
p2 = [p1 '/common'];
addpath (p2)
cd results/paramecium1/


% Importance Sampling Resampling I
disp('Importance Sampling Resampling with 1000 particles')

[m_k x_k res_ISR1] = IS_Resample(frames,bdryPts,ctrlPts,50,50,50,0.4,20,1000);

save res_ISR1

% Importance Sampling Resampling II
disp('Importance Sampling Resampling with 10000 particles')

[m_k x_k res_ISR2] = IS_Resample(frames,bdryPts,ctrlPts,50,50,50,0.4,20,10000);

save res_ISR2

% Importance Sampling - Resample Move
disp('Importance Sampling - Resample Move')

 [m_k x_k res_ISRM] = IS_ResampleMoveNtoM(frames,bdryPts,ctrlPts,50,50,50,0.4,20,1000,10,0.5,300,'RW');

save res_ISRM

% MAP
disp('MAP')

[m x res_MAP] = track_MAP(frames,bdryPts,ctrlPts,50,50,50,0.4,20,0.001,300,0);

save res_MAP


% save the results

disp('Saving the results')

% saving the raw frames

track_video_MCMC(frames,res_ISR1.output.m_k,res_ISR1.output.x_k,'ISR1','ISR1_')
disp('ISR1 saved')

track_video_MCMC(frames,res_ISR2.output.m_k,res_ISR2.output.x_k,'ISR2','ISR2_')
disp('ISR2 saved')

track_video_MCMC(frames,res_ISRM.output.m_k,res_ISRM.output.x_k,'ISRM','ISRM_')
disp('ISRM saved')

track_video(frames,res_MAP.output.m,res_MAP.output.x,'MAP','MAP_')
disp('MAP saved')

cd ..
cd ..
% removing the paths for the matlab functions
rmpath(p1)
rmpath(p2)
