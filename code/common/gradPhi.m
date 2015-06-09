function [G m z] = gradPhi(pts,ctrlPts,sig,dt,nofSteps,alpha);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
%
% This function calculates the gradient of Phi(.,1) with respect to alphas.
% Phi(.,tau) is the flow of a vector field with coefficients alpha.
% The function simultaneously constructs Phi(., tau) and it s gradient
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note that in this case bdry does not need to be only on the boundary.
   
   nofPts = size(pts,1);
   nofCtrlPts = size(ctrlPts,1);

   m = pts;
   z = ctrlPts;

   G = zeros(2*nofPts,2*nofCtrlPts); 

   for i=1:nofSteps
   
        % Building up the kernel for the boundary points
        dist_pts =  bsxfun(@minus,m(:,1),ctrlPts(:,1).').^2 + ...
                    bsxfun(@minus,m(:,2),ctrlPts(:,2).').^2;

        K_pts = exp(-dist_pts/(2*sig*sig));

        % Building up the kernel for the control points
	dist_ctrl = bsxfun(@minus,z(:,1),ctrlPts(:,1).').^2 + ...
                    bsxfun(@minus,z(:,2),ctrlPts(:,2).').^2;
        K_ctrl = exp(-dist_ctrl/(2*sig*sig));


	z2 = ctrlPts(:,1)*ones(1,nofPts);
	m2 = ones(nofCtrlPts,1)*(m(:,1))';
	dist1 =  z2 - m2 ;

	z2 = ctrlPts(:,2)*ones(1,nofPts);
	m2 = ones(nofCtrlPts,1)*(m(:,2))';
	dist2 =  z2 - m2 ;

        G11 = G(1:nofPts,1:nofCtrlPts);
        G21 = G(nofPts+1:end,1:nofCtrlPts);
        G12 = G(1:nofPts,nofCtrlPts+1:end);
        G22 = G(nofPts+1:end,nofCtrlPts+1:end);
 

	dist11 = (dist1'.*K_pts/(sig*sig))*alpha(:,1);
	dist22 = (dist2'.*K_pts/(sig*sig))*alpha(:,2);
	dist12 = (dist1'.*K_pts/(sig*sig))*alpha(:,2);
	dist21 = (dist2'.*K_pts/(sig*sig))*alpha(:,1);
       
        
        dist11 = dist11(:,ones(nofCtrlPts,1));
        dist22 = dist22(:,ones(nofCtrlPts,1));
        dist12 = dist12(:,ones(nofCtrlPts,1));
        dist21 = dist21(:,ones(nofCtrlPts,1));

        D1 = [dist11 dist11;dist12 dist12];
        D2 = [dist21 dist21;dist22 dist22];

        % evolving the points
        m(:,1) = m(:,1) + dt*K_pts*alpha(:,1);
        m(:,2) = m(:,2) + dt*K_pts*alpha(:,2);
        z(:,1) = z(:,1) + dt*K_ctrl*alpha(:,1); 
        z(:,2) = z(:,2) + dt*K_ctrl*alpha(:,2);

        G1_temp = [G11 G12;G11 G12];
        G2_temp = [G21 G22;G21 G22];

        % here I am adding G + blkdiag(K_pts,K_pts), it is faster since I am not allocating extra memory
         G(1:nofPts,1:nofCtrlPts) = G(1:nofPts,1:nofCtrlPts) + dt*K_pts;
         G(nofPts+1:end,1+nofCtrlPts:end) = G(1+nofPts:end,1+nofCtrlPts:end) + dt*K_pts;

        G = G + dt*(D1.*G1_temp + D2.*G2_temp); 

    end

