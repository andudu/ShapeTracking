function m_new = move_affine(m,par)
%calculate center of mass
%this could be moved outside so that I do not do it for each sample
%figure(2)
%plot(m(:,1),m(:,2))
%hold on
nofBdryPts = size(m,1);
cm_x = sum(m(:,1))/nofBdryPts;
cm_y = sum(m(:,2))/nofBdryPts;
m(:,1) = m(:,1) - cm_x;
m(:,2) = m(:,2) - cm_y;

  E_1 = [1 0 0;0 1 0; 0 0 0];
  E_2 = [1 0 0;0 -1 0; 0 0 0];
  E_3 = [0 -1 0;1 0 0; 0 0 0];
  E_4 = [0 1 0;1 0 0; 0 0 0];
  E_5 = [0 0 1;0 0 0;0 0 0];
  E_6 = [0 0 0; 0 0 1; 0 0 0];
  L = par(1)*E_1+par(2)*E_2+par(3)*E_3+par(4)*E_4+par(5)*E_5+par(6)*E_6; 


R = expm(L);
m_temp = [m ones(nofBdryPts,1)]';
m_new = (R*m_temp)';
m_new = m_new(:,1:2);
m(:,1) = m(:,1) + cm_x;
m(:,2) = m(:,2) + cm_y;
m_new(:,1) = m_new(:,1) + cm_x;
m_new(:,2) = m_new(:,2) + cm_y;
%plot(m_new(:,1),m_new(:,2))
