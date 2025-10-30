function loss=loss_fun_normal_Abeta_scalar(para, Age, Abeta, L, w)
y0=para(4)*ones(68, 1); % (68,1)
[t,y]=ode45(@(t,y) ODEmodel_original(para(1),para(2),para(3),L,y,t),Age,y0);  
err=y(ismember(t, Age(2:end-1)), :)' - Abeta;
loss = sum(sum(err.^2)./sum(Abeta.^2, 1)) + w*sum(sum((y(ismember(t,100),:)'-1).^2)); % relative mse

end