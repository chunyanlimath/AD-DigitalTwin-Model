function dydt=ODEmodel_original_68D(D,lambda_Abeta,K_Abeta,L,y,t) % the sign of each term in model is absorbed into param
% D>0 diffusion coeff, lambda>0 growth rate, K>0 carrying capacity
% D is a scalar, lambda and K are two 68D column vectors, L is graph Laplacian
% matrix with size 68 by 68.

dydt = -D*L*y+lambda_Abeta.*K_Abeta.*y-lambda_Abeta.*y.^2; % 

end
