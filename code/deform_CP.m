function [m_new x_i_new v]= deform_CP(m,x_i,x_i_new,alpha_k, Sigma, Sigma_0, dt)
nofBdryPts = size(m,1);
nofCtrlPts = size(x_i,1);

%Create the K_0 kernel matrix which moves the control points
z2 = x_i(:,1).^2*ones(1,nofCtrlPts);
m2 = ones(nofCtrlPts,1)*(x_i_new(:,1).^2)';
dist =  z2 - 2* x_i(:,1)*x_i_new(:,1)' + m2 ;
z2 = x_i(:,2).^2 * ones(1,nofCtrlPts) ;
m2 = ones(nofCtrlPts,1)*(x_i_new(:,2).^2)';
dist =  dist + z2 - 2* x_i(:,2)*x_i_new(:,2)' + m2 ;
K_0_ctrl = exp(-dist/(2*Sigma_0*Sigma_0))/(2*pi*Sigma_0*Sigma_0);

%Create the k_0 kernel matrix which moves the boundary points
z2 = x_i(:,1).^2*ones(1,nofBdryPts);
m2 = ones(nofCtrlPts,1)*(m(:,1).^2)';
dist =  z2 - 2* x_i(:,1)*m(:,1)' + m2;
z2 = x_i(:,2).^2 * ones(1,nofBdryPts);
m2 = ones(nofCtrlPts, 1) * (m(:,2).^2)';
dist =  dist + z2 - 2* x_i(:,2)*m(:,2)' + m2;
K_0_bdry = exp(-dist/(2*Sigma_0*Sigma_0))/(2*pi*Sigma_0*Sigma_0);

%Updating the control points
v(:,1) = K_0_ctrl'*alpha_k(:,1);
v(:,2) = K_0_ctrl'*alpha_k(:,2);
x_i_new(:,1) = x_i_new(:,1) + dt*v(:,1);
x_i_new(:,2) = x_i_new(:,2) + dt*v(:,2);

%Updating the boundary points
v = [];
v(:,1) = K_0_bdry'*alpha_k(:,1);
v(:,2) = K_0_bdry'*alpha_k(:,2);
m_new(:,1) = m(:,1) + dt*v(:,1);
m_new(:,2) = m(:,2) + dt*v(:,2);


function res = Kernel(x,y,Sigma)
    res = exp(-sum((x-y).^2)/(2*Sigma*Sigma))/(2*pi*Sigma*Sigma);
end

end