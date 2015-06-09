function tube = tubular_nbhd(image,bdryPts,width)
% I will scale it after that
nofBdryPts = size(bdryPts,1);
bdryPts = round(bdryPts);%when rounding I might start getting issues with boundaries
ind = sub2ind(size(image),bdryPts(:,2),bdryPts(:,1));
tube = zeros(size(image));
tube(ind)=1;
tube = bwmorph(tube,'dilate',width);
%  imagesc(tube)
%  hold on
%  plot(bdryPts(:,1),bdryPts(:,2),'Linewidth',2)