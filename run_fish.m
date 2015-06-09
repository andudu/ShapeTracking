 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fish sequence tracking.
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
disp('1 - fish sequence')
disp('2 - fish sequence with noise and grid')

reply = input('');
if (reply~=0&&reply~=1&&reply~=2) 
  error('Your input should be an integer between 0 and 4.')
end

a = input('Would you like to run the affine tracking? y/n ','s')

if (a~='y'&&a~='n')
  error('Your put should either y or n')
end


% loading the data
disp('Loading the data')
cd data
load 'fish_data.mat'
cd ..

ctrlPts = ctrlPts(2:end-1,:);

% creating a mat file to store the results
empty = [];
save 'results/fish/fish_mat' 'empty'


if (reply==0||reply==1)

% Experiment 1 - tracking the fish sequence 

disp('Tracking the fish sequence')

s = RandStream('mt19937ar','Seed',0)
RandStream.setGlobalStream(s)

if (a=='y')
  % affine tracking
  par_aff = [0.05 0 0.1 0 10 10];
  [m_average1 w_k P_k1] = track_affine_path2(frames(:,:,60:115),bdryPts,500,500,par_aff,{});
end
% nonlinear tracking
[m_k1 x_k1] = track_CP_sd(frames(:,:,60:115),bdryPts,ctrlPts,5000,5000,10,20,5,5,squeeze(mean(P_k1,3)),{});

% saving the matlab variables 
% save ('results/fish/fish_mat','res1','-append')
save ('results/fish/fish_mat','P_k1','m_k1','x_k1','-append')

% recording the results
track_video(frames(:,:,61:115),m_k1,x_k1,'results/fish/fish','results/fish/fish')


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Experiment 2 - Fish sequence + noise and occlusion

if (reply==0||reply==2)
disp('Tracking the fish sequence with noise and occlusion')


% resetting the random stream
s = RandStream('mt19937ar','Seed',0)
RandStream.setGlobalStream(s)

if (a == 'y')
  % affine tracking
  par_aff = [0.05 0 0.1 0 10 10];
  [m_average2 w_k P_k2] = track_affine_path2(frames_gridded(:,:,60:115),bdryPts,500,500,par_aff,{});
end
% nonlinear tracking
[m_k2 x_k2] = track_CP_sd(frames_gridded(:,:,60:115),bdryPts,ctrlPts,5000,5000,10,20,5,5,squeeze(mean(P_k2,3)),{});


% saving the matlab variables 
% save ('results/fish/fish_mat','res2','-append')
save ('results/fish/fish_mat','P_k2','m_k2','x_k2','-append')

% recording the results
track_video(frames_gridded(:,:,61:115),m_k2,x_k2,'results/fish/fish_gridded','results/fish/fish_gridded')

end



rmpath(p1)
rmpath(p2)
matlabpool close
