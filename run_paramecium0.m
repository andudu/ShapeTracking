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
load('param.mat')

cd ..
% run the experiments 
p = pwd;
p1 = [p '/code'];
p2 = [p1 '/common'];
addpath (p1)
addpath (p2)

% Affine Tracking
par_aff = [0.01 0.01 0.05 0.01 2 2];
% m_k_affine = IS_Resample_affine(frames,bdryPts,5000,par_aff);
m_affine = track_affine_path1(frames,bdryPts,1000,1000,par_aff,[]);


disp('Saving the results')

track_video(frames,m_affine,[],'results/paramecium0/param_affine','results/paramecium0/param_affine')
% track_video(frames,m_k_affine,[],{},{})

save ('results/paramecium0/m_affine', 'm_affine');


disp('ISR Affine saved')



% Importance Sampling Resampling I
disp('Importance Sampling Resampling with 1000 particles')
[m_k x_k res_ISR] = IS_Resample(frames,bdryPts,ctrlPts,50,50,50,0.1,15,1000);
save ('results/paramecium0/res_ISR','res_ISR'); 

% save the results
disp('Saving the results')
track_video(frames,squeeze(mean(m_k,3)),squeeze(mean(x_k,3)),'results/paramecium0/param_deform','results/paramecium0/param_deform')
disp('ISR saved')


rmpath(p1)
rmpath(p2)
matlabpool close