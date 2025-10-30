function dydt=ODEmodel_original(D,lambda_Abeta,K_Abeta,L,y,t) % the sign of each term in model is absorbed into param
% D>0 diffusion coeff, lambda>0 growth rate, K>0 carrying capacity
% D, lambda, K are 3 scalars and L is graph Laplacian matrix with shape 68
% by 68 and y is a 68D vector.  
dydt = -D*L*y+lambda_Abeta*K_Abeta*y-lambda_Abeta*y.^2;
end
