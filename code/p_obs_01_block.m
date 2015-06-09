function p = p_obs_01_block(V0,m0,par_obs,type)

%this will cause an issue if the object is close to the boundary (need to fix it)
a1 = max(round(par_obs(11,1) - par_obs(8,1)/2),1);
a2 = max(round(par_obs(10,1) - par_obs(9,1)/2),1);
b1 = max(round(par_obs(11,1)) + round(par_obs(8,1)/2),size(V0,1));
b2 = max(round(par_obs(10,1)) + round(par_obs(9,1)/2),size(V0,2));


V0 = V0(a1:b1,a2:b2);
m0(:,1) = m0(:,1) - a2;
m0(:,2) = m0(:,2) - a1;

% check if the curve is inside the region
% if not project the outlying values inside

%  m0(find(m0<1)) = 1;
%  m0(find(m0(:,1)>size(V0,2)),1) = size(V0,2);
%  m0(find(m0(:,2)>size(V0,1)),2) = size(V0,1);

% check if the curve is inside the region
% if not remove the outlying points

m0(find(m0(:,1)<1),:) = [];
m0(find(m0(:,2)<1),:) = [,];
m0(find(m0(:,1)>size(V0,2)),:) = [];
m0(find(m0(:,2)>size(V0,1)),:) = [];


%  imagesc(V0)
%  hold on
%  plot(m0(:,1),m0(:,2))
%  plot(par_obs(10,1),par_obs(11,1),'*');


L = roipoly(V0,m0(:,1),m0(:,2));
%L = scanPolyCurve(m0, [size(V0,2), size(V0,1)])' ;
%L(L==0) = 1;
% L = L(1:2:end,1:2:end);
% V0 = V0(1:2:end,1:2:end);

if strcmp(type,'gaussian')
  cin = par_obs(1,1); cout = par_obs(2,1);
  varin = par_obs(3,1); varout = par_obs(4,1);
  kin = par_obs(5,1); kout = par_obs(6,1);
  p_max = par_obs(7,1);

          kin=1;
        kout=1;


  n = size(V0,1)*size(V0,2);
  exp_1 = -sum((V0(L==1) - cin).^2)/(2*varin)-sum(sum(L==1))*log(2*pi*varin)/2;
  exp_2 = -sum((V0(L==0) - cout).^2)/(2*varout)-sum(sum(L==0))*log(2*pi*varout)/2;
  % exp_1 = -sum((V0(L==1) - cin).^2)/(2*varin);
  % exp_2 = -sum((V0(L==0) - cout).^2)/(2*varout);
%  
%    disp(log(kin))
%    disp(log(kout))
 % exp_1 = -sum((V0(L==1) - cin).^2)/(2*varin)-sum(sum(L==1))*log(kin);
 % exp_2 = -sum((V0(L==0) - cout).^2)/(2*varout)-sum(sum(L==0))*log(kout);
%    disp(exp_1)
%    disp(exp_2)
 %  p_max=0;
  p = exp_1+exp_2 + n*log((kout+kout)/2);
 % p = exp(exp_1+exp_2 +p_max);
elseif strcmp(type,'beta')

  cin = par_obs(1,1); cout = par_obs(2,1);
  varin = par_obs(3,1); varout = par_obs(4,1);
  alphain = (1-cin)*cin*cin/varin - cin;
  alphaout = (1-cout)*cout*cout/varout - varout;
  betain = alphain*(1-cin)/cin;
  betaout = alphaout*(1-cout)/cout;
  disp(alphain)
  disp(betain)

  in = sum(sum(L==2));
  out = sum(sum(L==1));
  p = (prod(V0(L==2))*prod(V0(L==1)))^(alphain-1)*(prod(1-V0(L==2))*prod(1-V0(L==1)))^(betain-1)/(beta(alphain,betain)^in*beta(alphaout,betaout)^out);
end
