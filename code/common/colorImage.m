function outImage = colorImage(inImage,bdryPts,sizes,mu,sig)
  outBdryPts = bdryPts(1:sizes(1,1),:);
  inBdryPts = bdryPts(sizes(1,1)+1:end,:);
  BW= roipoly(inImage,outBdryPts(:,1),outBdryPts(:,2));
  outImage = mu(2,1)*BW+(1-BW)*mu(1,1);
  BW = roipoly(zeros(size(inImage)), inBdryPts(:,1),inBdryPts(:,2));
  outImage = (1-BW).*outImage + BW*mu(3,1);
  outImage = outImage + sig*randn(size(inImage));
end