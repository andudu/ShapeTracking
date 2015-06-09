function logL = logObsL_region2(V0,m0,sizes,par_obs,type)
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
  a1 = max(round(par_obs(11,1) - par_obs(10,1)/2),1);
  a2 = max(round(par_obs(10,1) - par_obs(9,1)/2),1);
  b1 = min(round(par_obs(11,1)) + round(par_obs(8,1)/2),size(V0,1));
  b2 = min(round(par_obs(10,1)) + round(par_obs(9,1)/2),size(V0,2));

  % new image + new curve
  V0 = V0(a1:1:b1,a2:1:b2);
  m0(:,1) = m0(:,1) - a2;
  m0(:,2) = m0(:,2) - a1;

end

% check if the curve is inside the region
% if not remove the outlying points
% (this part is buggy, it might happen that the whole curve is outside so there are no pts left)
m0(find(m0(:,1)<1),:) = [];
m0(find(m0(:,2)<1),:) = [,];
m0(find(m0(:,1)>size(V0,2)),:) = [];
m0(find(m0(:,2)>size(V0,1)),:) = [];



bdryPts1 = m0(1:sizes(1,1),:);
bdryPts2 = m0(sizes(1,1)+1:end,:);

L = colorImage(V0,m0,sizes,[1 2 3]',0);


L = L(1:1:end,1:1:end);
V0 = V0(1:1:end,1:1:end); % if I scale my images I need to scale p_max too

%    mu(1,1) = mean(V0(L==1));
%    mu(2,1) = mean(V0(L==2));
%    mu(3,1) = mean(V0(L==3));

      mu(1,1) = par_obs(1,1);
      mu(2,1) = par_obs(2,1);
      mu(3,1) = par_obs(3,1);

  var = par_obs(4:6,1);
  p_max = par_obs(7,1);

  n = size(V0,1)*size(V0,2);
%    exp_1 = -sum((V0(L==1) - cin).^2)/(2*varin)-sum(sum(L==1))*log(2*pi*varin)/2;
%    exp_2 = -sum((V0(L==0) - cout).^2)/(2*varout)-sum(sum(L==0))*log(2*pi*varout)/2;
  E_1 = -sum((V0(L==1) - mu(1,1)).^2)/(2*var(1,1));
  E_2 = -sum((V0(L==2) - mu(2,1)).^2)/(2*var(2,1));
  E_3 = -sum((V0(L==3) - mu(3,1)).^2)/(2*var(3,1));

  E_4 = -sum(sum(L==1))*log(sqrt(var(1,1))) - sum(sum(L==2))*log(sqrt(var(2,1)))-sum(sum(L==3))*log(sqrt(var(3,1)));
%    E_4 = 0;
%    exp_1 = -sum((V0(L==1) - cin).^2)/(2*varin)-sum(sum(L==1))*log(kin);
%    exp_2 = -sum((V0(L==0) - cout).^2)/(2*varout)-sum(sum(L==0))*log(kout);
%    disp(exp_1+exp_2)
 %  p_max=0;
%  n
%  
%      disp(exp_1+exp_2 +exp_3)
%    p = exp(exp_1+exp_2 + exp_3);
%     p = exp(exp_1+exp_2 +exp_3-0.3*p_max);
logL = E_1 + E_2 + E_3 + E_4;