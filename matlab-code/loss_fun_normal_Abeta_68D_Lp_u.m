function loss=loss_fun_normal_Abeta_68D_Lp_u(uv_Lp, para, u0_para, Age,Abeta,Lp, w, w_u, w_su, w_v, w_sv)
tspan=Age;
u_Lp = uv_Lp(1:68);  % （68， 1）
v_Lp = uv_Lp(69:136); %（68， 1）
L = Lp + 0.5*[ones(1,68)*v_Lp * diag(u_Lp) + ones(1, 68)*u_Lp*diag(v_Lp)] - 0.5*(u_Lp*v_Lp'+v_Lp*u_Lp');
% Save L into L.mat
% save('L.mat', 'L');
% options = odeset('RelTol',1e-8, 'AbsTol',1e-10);
[t,y]=ode45(@(t,y) ODEmodel_original_68D(para(1),para(2:69),para(70:137),L,y,t),tspan,u0_para);  

err=y(ismember(t, tspan(2:end-1)), :)' - Abeta;
loss = sum(sum(err.^2)./sum(Abeta.^2, 1)) + w*sum(sum((y(ismember(t,100),:)'-1).^2)) + w_u * (norm(u_Lp)-1)^2 + w_su * sum(abs(u_Lp)) + w_v*(norm(v_Lp)-1)^2 +w_sv*sum(abs(v_Lp));
end

% ODE Output Function
% function status = odeOutput(t, y, flag)
%     persistent dydtLog;
%     if isempty(dydtLog), dydtLog = {}; end
% 
%     if strcmp(flag, 'init')
%         dydtLog = {};
%     elseif strcmp(flag, '')
%         % Store intermediate values
%         dydtLog{end+1} = struct('t', t, 'y', y);
%     elseif strcmp(flag, 'done')
%         % Save to a file or analyze
%         save('debug_dydtLog.mat', 'dydtLog');
%     end
%     status = 0;
% end