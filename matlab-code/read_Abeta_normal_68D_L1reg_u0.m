% fix parameters obtained from 68D+L1 reg (eps) selected best based on fit loss
% and optimize u0 (2nd level optimization) to improve the prediction performance by minimizing fit loss only.
% note that the shape of solution should not changed since the para is fixed!
% Note that the discrete population Laplacian is used as the L_p. 

% load Lp matrix with size 68 * 68. 
load('Lp68.mat');
% ------- Load the required data --------------------------------------------------------------------
data = readmatrix('With_age_UCBERKELEYAV45_04_26_22.csv');
Rid = data(:, 1);
PID = unique(Rid);
Age = data(:, 4);

Abeta_ref = data(:, 12);
Abeta = data(:, 42:2:181)./ Abeta_ref;
columns_with_all_nan = all(isnan(Abeta), 1); % Calculate the indices of columns that are entirely NaN
Abeta = Abeta(:, ~columns_with_all_nan); % Remove those columns
disp('Size of Abeta after removing columns with all NaN values:');
disp(size(Abeta)); % (3086, 68)

% --------- normalization ----------------
Abeta_Max = max(max(Abeta))
Abeta = Abeta / Abeta_Max;

best_fit = [];
Num_pts = [];

% -------- initial setting ---------------------------------------------------------------------------
plots_folder = 'plots_Abeta_68D/fmincon_Lp_68D_regL1_vary_wreg/';
data_filename = fullfile(plots_folder, 'collect_best_w_reg_para_fitloss_Sshape.xlsx');
para0_all = readmatrix(data_filename); % (471, 137)

% Generate the folder and filename for saving the results
plots_folder = 'plots_Abeta_68D_u0/fmincon_nolsq_w_Lp_68D_regL1_u0/';   
if ~exist(plots_folder, 'dir')
    mkdir(plots_folder);
end  

final_1d=readmatrix('plots_Abeta_scalar/fmincon_nolsq_w_0_00015625Lp/summary_best_w_patient_record.xlsx');
pid_idx = final_1d(:, 1);
w_all = final_1d(:, 2);
Abeta0_all = final_1d(:, 6);

parfor i = 1:length(pid_idx)
    index = find(Rid == Rid(pid_idx(i)));
    % Define log file for the current patient
    log_filename = fullfile(plots_folder, ['patient_log_' num2str(index(end)) '.txt']);
    fid = fopen(log_filename, 'w'); % Open the log file in write mode
    
    try
        if length(index) > 2 
            patient_record = []; % param for each patient for each random initial guess
            patient_grad = [];
            u0_all = [];
            disp(['====== PID is ' num2str(index(end))]);
            L = Lp;
            w = w_all(i);
            para = para0_all(i, :)'; % (137, 1) fixed 137 para (theta) obtained from 68D+L1 reg(eps)
            tspan = [50; Age(index(1:end-1)); 100];
            tspan1 = [50; Age(index(1:end)); 100];
            u0 = Abeta0_all(i)*ones(68, 1); % (68, 1)
            u0_initial = zeros(68, 1);  % eps0=u0_real-u0
            lb = -u0;   % (68, 1) lower bound for u0=0
            ub = 1 * ones(68, 1);  % upper bound for u0=1
            w1 = 1e3;
            for k=1:11
                disp(['w1=', num2str(w1)]);
                % ------- optimize u0 -------
                [u0_para, fval, exitflag, output, lambda_out, grad, ~] = fmincon(@(u0_para) loss_fun_normal_Abeta_68D_regL1_u0(u0_para, u0, para, tspan, Abeta(index(1:end-1), :)', double(L), w, w1), u0_initial, [], [], [], [], lb, ub, [], optimset('MaxFunEvals', 1e7));

                grad_vec = grad(:);
                patient_grad = [patient_grad, grad_vec]; % (68, 471)
                disp(['exitflag=', num2str(exitflag)]);

                [t, y] = ode45(@(t, y) ODEmodel_original_68D(para(1), para(2:69), para(70:137), double(L), y, t), tspan1 , u0_para+u0);

                err = y(ismember(t, Age(index(1:end))), :)' - Abeta(index(1:end), :)';
                pred_err_lsq = sqrt(sum(err.^2, 1)) ./ sqrt(sum(Abeta(index(1:end), :).^2, 2))'; % sum over rows, keep column
                fit_loss = sqrt(sum(pred_err_lsq(1:end-1).^2)/(length(index)-1));
           
                err_T = sqrt(sum((y(ismember(t, 100), :)' - 1).^2));

                patient_record = [patient_record; [index(end); w1; fval; fit_loss; pred_err_lsq(end); exitflag; pred_err_lsq'; err_T; u0_para+u0]'];
                u0_all=[u0_all, u0_para+u0];
                disp(['fmincon err is ', num2str(pred_err_lsq)]);
                disp(['fmincon err for T=100 is ', num2str(err_T)]);
                if w1>1.0
                    w1 = w1/10;
                else
                    w1=w1/2;
                end
            end
            w_str = sprintf('%.8f', w1); 
            w_str = strrep(w_str, '.', '_');

            data_filename = fullfile(plots_folder, ['fmincon_grad_PID' num2str(index(end)) '_w' w_str '.csv']);
            writematrix(patient_grad, data_filename);

            data_filename = fullfile(plots_folder, ['patient_record' num2str(index(end)) '_w' w_str '.csv']);
            writematrix(patient_record, data_filename);

            data_filename = fullfile(plots_folder, ['patient_record_u0' num2str(index(end)) '_w' w_str '.csv']);
            writematrix(u0_all, data_filename);
        end
    catch ME
        % Log error information
        fprintf(fid, 'Error occurred: %s\n', ME.message);
    end
    
    % Close the log file
    fclose(fid);        
end





