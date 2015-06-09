function logL = logObsL_region(V0,m0,par_obs,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% E = logObsL_region(V0,mo,par_obs,type)
%
% logObsL_region - is a function calclulating the logarithm of the observation likelihood based on 
% 			a region based model of the object in the image: each pixel is assumed 
%  			a normally distributed random variable  with mean and variance depending 
%			on being inside or outside of the contour
%	
% Inputs: V0 - the image, i.e. the observation
% 	  m0 - a nofBdryPts times 2 array which contains the coordinates of the boundary
%         par_obs - all the parameters needed to describe the normal densities
%         
% OUTPUT: loL - a double value containing the log of the probability
%         
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(type,'block')
  % calculating the new image dimensions
  % this will cause an issue if the object is close to the boundary (need to fix it)
  a1 = max(round(par_obs(9,1) - par_obs(6,1)/2),1);
  a2 = max(round(par_obs(8,1) - par_obs(7,1)/2),1);
  b1 = min(round(par_obs(9,1)) + round(par_obs(6,1)/2),size(V0,1));
  b2 = min(round(par_obs(8,1)) + round(par_obs(7,1)/2),size(V0,2));

  % new image + new curve
  V0 = V0(a1:1:b1,a2:1:b2);
  m0(:,1) = m0(:,1) - a2;
  m0(:,2) = m0(:,2) - a1;
end



% check if the curve is inside the region
% if not remove the outlying points

m0(find(m0(:,1)<1),:) = [];
m0(find(m0(:,2)<1),:) = [,];
m0(find(m0(:,1)>size(V0,2)),:) = [];
m0(find(m0(:,2)>size(V0,1)),:) = [];


% Label the image:  1 inside the contour; 0 outside the contour.
L = roipoly(V0,m0(:,1),m0(:,2));

% I can use the following lines to construct a sparser image
L = L(1:1:end,1:1:end);
V0 = V0(1:1:end,1:1:end);

% setting likelihood parameters
cin = par_obs(1,1); cout = par_obs(2,1);
varin = par_obs(3,1); varout = par_obs(4,1);

% Calculating the energy inside and outside
E_1 = sum((V0(L==1) - cin).^2)/(2*varin);
E_2 = sum((V0(L==0) - cout).^2)/(2*varout);
E_3 = sum(sum(L==1))*log(sqrt(varin))+sum(sum(L==0))*log(sqrt(varout));

logL = - E_1 - E_2 - E_3;
