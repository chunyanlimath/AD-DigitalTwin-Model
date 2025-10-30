function loss=loss_fun_normal_Abeta_68D_regL1_u0(u0_para, u0, para, Age, Abeta, L, w, w1)
    y0=u0_para+u0;
    [t,y]=ode45(@(t,y) ODEmodel_original_68D(para(1),para(2:69),para(70:137),L,y,t),Age,y0);  % u0_para:（68,1）
    err=y(ismember(t, Age(2:end-1)), :)' - Abeta;
    loss = sum(sum(err.^2)./sum(Abeta.^2, 1)) + w*sum(sum((y(ismember(t,100),:)'-1).^2))+w1*sum(abs(u0_para));
     
end  