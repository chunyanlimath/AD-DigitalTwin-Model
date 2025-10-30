function loss=loss_fun_normal_Abeta_68D_regL1(epsilon, theta0, Age,Abeta, Abeta0, L, w, w_reg)
    para=theta0+epsilon;
    [t,y]=ode45(@(t,y) ODEmodel_original_68D(para(1),para(2:69),para(70:137),L,y,t),Age, Abeta0);  
    err=y(ismember(t, Age(2:end-1)), :)' - Abeta; % y size=(num_t, 68)
    loss = sum(sum(err.^2)./sum(Abeta.^2, 1)) + w*sum(sum((y(ismember(t,100),:)'-1).^2)) + w_reg * sum(abs(epsilon));
     
end