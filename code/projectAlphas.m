function alpha_k_projected = projectAlphas(x_average,x_k,alpha_k,sig)

alpha_k_projected = zeros(size(alpha_k));
for k=1:size(x_k,3)
  K10 = K_matrix(x_average,x_k(:,:,k),sig^2);
  K11 = K_matrix(x_average,x_average,sig^2);
  alpha_k_projected(:,:,k) = inv(K11)*K10*alpha_k(:,:,k);
end
