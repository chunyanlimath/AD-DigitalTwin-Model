% in this code, we will implement regularization to overcome the overfitting issue for 68D case. 
% set theta=theta+epslion and we impose a regularization term $||\epslion||_2^2$ in loss function to minimize
% the initial of para (D, lambda, K) are given by the scalar baseline model
% D=D_scalar * ones(68, 1), etc.
% L = L_p, and L_p is the mean of all Laplacian matrix 
% Abeta0 = u0*ones(68, 1) where u0 is obtained from scalar case 

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

% -------- Normalization ----------------------------------------------------------------------------
Max_Abeta = max(max(Abeta));
Abeta = Abeta / Max_Abeta;

plots_folder = 'plots_Abeta_68D/fmincon_Lp_68D_regL1_vary_wreg/';
if ~exist(plots_folder, 'dir')
    mkdir(plots_folder);
end
% -------- initial setting (result of scalar baseline model obtained from 67 patients)---------------------------------------------------------------------------
final_1d=readmatrix('plots_Abeta_scalar/fmincon_nolsq_w_0_00015625Lp/summary_best_w_patient_record.xlsx');
best_fit = [];
Num_pts = [];
parfor i = 1:size(final_1d, 1)
    w = final_1d(i, 2);
    D = final_1d(i, 3);  % D = D * ones(68, 1)
    lambda = final_1d(i, 4);
    K = final_1d(i, 5);
    Abeta0 = final_1d(i, 6)*ones(68,1);
    index = find(Rid == Rid(final_1d(i, 1)));
    
    % Define log file for the current patient
    log_filename = fullfile(plots_folder, ['patient_log_' num2str(index(end)) '.txt']);
    fid = fopen(log_filename, 'w'); % Open the log file in write mode
    
    try
        if length(index) > 2
            patient_record = []; % param for each patient for each random initial guess
            patient_grad = [];
            epsilon_all = [];
            disp(['====== PID is ' num2str(index(end))]);

            L = Lp;

            lb = [-D; -lambda * ones(68, 1); 1-K * ones(68, 1)];  % (137, 1)
            ub = [1.0; 1 * ones(68, 1); 1.0 * ones(68, 1)]; % ub for epsilon in [0, 1]

            theta0 = [D; lambda * ones(68, 1); K * ones(68, 1)]; % (137, 1)
            epsilon0 = zeros(137, 1);
            tspan = [50; Age(index(1:end-1)); 100]; % (num_t, 1)
            tspan1 = [50; Age(index(1:end)); 100]; % (num_t+1, 1)
            w_reg = 1e3;
            for k=1:11
                disp(['w_reg=', num2str(w_reg)]);
                [epsilon, fval, exitflag, output, lambda_out, grad, ~] = fmincon(@(epsilon) loss_fun_normal_Abeta_68D_regL1(epsilon, theta0, tspan, Abeta(index(1:end-1), :)', Abeta0, double(L), w, w_reg), epsilon0, [], [], [], [], lb, ub, [], optimset('MaxFunEvals', 1e7));

                epsilon = double(epsilon);
                grad_vec = grad(:);
                patient_grad = [patient_grad, grad_vec];

                disp(['exitflag=', num2str(exitflag)]);
                para = epsilon + theta0;
                [t, y] = ode45(@(t, y) ODEmodel_original_68D(para(1), para(2:69), para(70:137), double(L), y, t), tspan1 , Abeta0);
                % y size (num_t, 68)
                err = y(ismember(t, Age(index(1:end))), :)' - Abeta(index(1:end), :)';
                pred_err_lsq = sqrt(sum(err.^2, 1)) ./ sqrt(sum(Abeta(index(1:end), :).^2, 2))'; % sum over rows, keep column
                fit_loss = sqrt(sum(pred_err_lsq(1:end-1).^2)/(length(index)-1)); 

                err_T = sqrt(sum((y(ismember(t, 100), :)' - 1).^2));
                patient_record = [patient_record; [index(end); w_reg; fval; fit_loss; pred_err_lsq(end); exitflag; pred_err_lsq'; err_T; para]'];
                epsilon_all=[epsilon_all, epsilon];
                disp(['fmincon err is ', num2str(pred_err_lsq)]);
                disp(['fmincon err for T=100 is ', num2str(err_T)]);
                % Reduce w by a factor of 1/2 and update para0
                if w_reg>1.0
                    w_reg = w_reg/10;
                else
                    w_reg=w_reg/2;
                end
                epsilon0 = epsilon; % homotopy idea for w
            end
            w_str = sprintf('%.8f', w_reg); 
            w_str = strrep(w_str, '.', '_');

            data_filename = fullfile(plots_folder, ['fmincon_grad_PID' num2str(index(end)) '_w' w_str '.csv']);
            writematrix(patient_grad, data_filename);

            data_filename = fullfile(plots_folder, ['patient_record' num2str(index(end)) '_w' w_str '.csv']);
            writematrix(patient_record, data_filename);

            data_filename = fullfile(plots_folder, ['patient_record_epsilon' num2str(index(end)) '_w' w_str '.csv']);
            writematrix(epsilon_all, data_filename);
        end
    catch ME
        % Log error information
        fprintf(fid, 'Error occurred: %s\n', ME.message);
    end
    
    % Close the log file
    fclose(fid);
end




