function par_obs = parObsL2(V0,bdryPts,sizes,W,H,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% par_obs = parObsL2(V0,bdryPts,sizes,W,H,type)
%
% This function initializes the parameters of the observation likelihood based on 
% the initial frame (it assumes the initial boundary correctly segments the region 
% of interest).
% 
% INPUT:
% 	V0 - a matrix representing the first frame in the video sequence
% 	bdryPts - a (nofBdryPts,2) arrray containing the locations of the 
%		boundary points of the object in the first frame
%       W - an integer specifying the width of the block for a block observation likelihood
%	H - an integer specifying the height of the block for a block observation likelihood
%
% OUTPUT:
%	par_obs - a variable containing all the parameters of needed to compute the 
%		observation likelihood
%       par_obs(1,1) = cin; par_obs(2,1) = cout; par_obs(3,1) = 
%	par_obs(4,1) = varin; par_obs(5,1) = varout; par_obs(6,1) =  
%	par_obs(8,1) = p_max;
%	par_obs(8,1) = H; par_obs(,1) = W;         
%       par_obs(10,1) - the x-axis coordinate for the center of mass
%	par_obs(11,1) - the y-axis coordinate for the center of mass
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


par_obs(8,1) = H; par_obs(9,1) = W; %size of the block (need to check the order)
% center of mass of the initial segmentation (in the later stages will use the center of mass of the previous mean curve) 
par_obs(10,1) = mean(bdryPts(:,1)); par_obs(11,1) = mean(bdryPts(:,2)); 

if strcmp(type,'block')
  % calculating the new image dimensions
  % this will cause an issue if the object is close to the boundary (need to fix it)
  a1 = max(round(par_obs(11,1) - par_obs(8,1)/2),1);
  a2 = max(round(par_obs(10,1) - par_obs(9,1)/2),1);
  b1 = min(round(par_obs(11,1)) + round(par_obs(8,1)/2),size(V0,1));
  b2 = min(round(par_obs(10,1)) + round(par_obs(9,1)/2),size(V0,2));

  % new image + new curve
  V0 = V0(a1:1:b1,a2:1:b2);
  bdryPts(:,1) = bdryPts(:,1) - a2;
  bdryPts(:,2) = bdryPts(:,2) - a1;
end


L = colorImage(V0,bdryPts,sizes,[1 2 3]',0);




par_obs(1,1) = mean(V0(L==1));
par_obs(2,1) = mean(V0(L==2));
par_obs(3,1) = mean(V0(L==3));

par_obs(4,1) = var(V0(L==1));
par_obs(5,1) = var(V0(L==2));
par_obs(6,1) = var(V0(L==3));


%par_obs(4,1) = 0.1;
%par_obs(5,1) = 1;
%par_obs(6,1) = 1;

%par_obs(7,1) = -sum(sum((L.*(V0 - cin)).^2/(2*varin) + ((1-L).*(V0-cout)).^2/(2*varout)));
%par_obs(3,1) = 0.1; par_obs(4,1) = 0.1;
