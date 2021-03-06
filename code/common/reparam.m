function  [bdryPts, ctrlPts] = reparam(bdryPts,ctrlPts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   function  [bdryPts, ctrlPts] = reparam(bdryPts,ctrlPts)
%  
%   This function reparameterizes the curve bdryPts according to
%   arclength. It also reparameterizes the control points and 
%   ensures that the number of control points stays the same.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin<2)
  ctrlPts = [];
end

    nofBdryPts = size(bdryPts,1);
    nofCtrlPts = size(ctrlPts,1);
    % reparameterizing
    bdryPts = arcLength(bdryPts,1/nofBdryPts);
    if (size(bdryPts,1)~= nofBdryPts)
      disp('Error in curve reparameterization.')
    end

    if (~isempty(ctrlPts))
      % maybe it should be modified so that I don't end up with the first and last control points being equal
      ctrlPts = bdryPts(1:(floor((nofBdryPts-1)/(nofCtrlPts-1))):end,:);
      if (size(ctrlPts,1)<nofCtrlPts)
	disp('Error in control points reparameterization')
      end

      if (size(ctrlPts,1)~=nofCtrlPts)
	ctrlPts = ctrlPts(1:nofCtrlPts,:);
      end
    end
end
