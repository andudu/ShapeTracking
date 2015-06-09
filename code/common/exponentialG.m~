function [m z m_path z_path] = exponentialG(m0,z0,v0,sig,dt,nofSteps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% [m z m_path z_path] = exponentialG(m0,z0,v0,sig,dt,nofSteps)
% 
% This function evolves a set of boundary points m0 and a set of cotrol points z0
% along the group exponential map with initial momentum v0. The width of the kernel is sig,
% dt is the size of the time step and nofSteps is the number of steps in the evolution.
% 
% 
% OUTPUT:
%       m - the final position of the curve
%       z - the final position of the control points
%       m_path - stores the curve at each intermediate step (nofBdryPts,2,nofSteps)
%       z_path - stores the ctrlPts at each intermediate step (nofCtrlPts,2,nofSteps)
%
% Notes: I should make a switch to save the paths only if needed
%        I should make a switch to plot the curves only if needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    m = m0;
    z = z0;
    m_path = zeros(size(m0,1),2,nofSteps);
    z_path = zeros(size(z0,1),2,nofSteps);
%      plot(m(:,1),m(:,2),'r','Linewidth',2)
    for i=1:nofSteps
        % Building up the kernel for the boundary points
	dist = bsxfun(@minus,m(:,1),z0(:,1).').^2 + ...
                    bsxfun(@minus,m(:,2),z0(:,2).').^2;
	K = exp(-dist/(2*sig*sig));

        % Building up the kernel for the control points
	dist = bsxfun(@minus,z(:,1),z0(:,1).').^2 + ...
                    bsxfun(@minus,z(:,2),z0(:,2).').^2;
	K00 = exp(-dist/(2*sig*sig));



        % evolving the points
        m(:,1) = m(:,1) + dt*K*v0(:,1);
        m(:,2) = m(:,2) + dt*K*v0(:,2);
        z(:,1) = z(:,1) + dt*K00*v0(:,1);
        z(:,2) = z(:,2) + dt*K00*v0(:,2);
        m_path(:,:,i) = m; 
        z_path(:,:,i) = z;
%          plot(m(:,1),m(:,2),'b.','Linewidth',0.1)
%          pause(0.01)
    end
%      plot(m(:,1),m(:,2),'r','Linewidth',2)
end