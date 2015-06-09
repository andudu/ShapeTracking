function [m_new x_new] = move_CP_sd(m,x,P)
%calculate center of mass
%this could be moved outside so that I do not do it for each sample
%figure(2)
%plot(m(:,1),m(:,2))
%hold on
nofBdryPts = size(m,1);
nofCtrlPts = size(x,1);

m_temp = [m ones(nofBdryPts,1)]';
m_new = (P*m_temp)';
m_new = m_new(:,1:2);
x_temp = [x ones(nofCtrlPts,1)]';
x_new = (P*x_temp)';
x_new = x_new(:,1:2);

%plot(m_new(:,1),m_new(:,2))
%calculate center of mass
%this could be moved outside so that I do not do it for each sample
%figure(2)
%plot(m_new(:,1),m_new(:,2))
%hold on
