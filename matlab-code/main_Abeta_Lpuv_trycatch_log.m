% load Lp matrix with size 68 * 68. 
load('Lp68.mat');
load('A.mat');
Lp_const = parallel.pool.Constant(Lp);
% Broadcast A to all workers using parallel.pool.Constant
A_const = parallel.pool.Constant(A);

% ---------- load data ------------------ 
% load the Abeta data in all brain regions for all patients collected
data=readmatrix('With_age_UCBERKELEYAV45_04_26_22.csv');
Rid=data(:,1); % research participant identifier which is unique for every patient
PID=unique(Rid); % one patient might have multiple fMRI scan or PET scan 
Age=data(:,4);
% normalization 
Abeta_ref = data(:, 12); % COMPOSITE_REF_SUVR 
Abeta=data(:,42:2:181)./Abeta_ref;  % Abeta info in Cortext region
% remove two unknown columns 
columns_with_all_nan = all(isnan(Abeta), 1); % Calculate the indices of columns that are entirely NaN
Abeta = Abeta(:, ~columns_with_all_nan); % Remove those columns
disp('Size of Abeta after removing columns with all NaN values:');
disp(size(Abeta)); % (3086, 68)

% --------- normalization ----------------
Abeta_Max = max(max(Abeta))
Abeta = Abeta / Abeta_Max;

plots_folder = 'plots_Abeta_68D/fmincon_Lp_68D_regL1_vary_wreg/';
data_filename = fullfile(plots_folder, 'collect_best_w_reg_para_fitloss_Sshape.xlsx');
para0_all=readmatrix(data_filename); % (471, 137)

% Generate the folder and filename for saving the results
plots_folder = 'plots_Abeta_68D_u0/fmincon_nolsq_w_Lp_68D_regL1_u0/';  %_homotopy_w
data_filename = fullfile(plots_folder, 'summary_patient_record_u0.xlsx');
u0_all = readmatrix(data_filename);  % (471, 68)

final_1d=readmatrix('plots_Abeta_scalar/fmincon_nolsq_w_0_00015625Lp/summary_best_w_patient_record.xlsx');
pid_idx = final_1d(:, 1);
w_all = final_1d(:, 2);
% Generate the folder and filename for saving the results
plots_folder = 'plots_Abeta_Lp_uu/fmincon_nolsq_68D_u0_Lp_lowrank_u_wu1000_ws1000-cluster/';  %_homotopy_w
if ~exist(plots_folder, 'dir')
    mkdir(plots_folder);
end  

w_u = 10000;  % weight for unit norm of u.
w_v = 10000;  % weight for unit norm of v.
Accuracy=[];
best_fit=[];
Num_pts=[];  % 

A(A == 0) = NaN;
% Find the minimum for each row ignoring NaN
minNonzeroValues = min(A, [], 2, 'omitnan');
maxNonzeroValues = max(A, [], 2, 'omitnan');
ub = double([1-minNonzeroValues; 1-minNonzeroValues]);
lb = double([-maxNonzeroValues; -maxNonzeroValues]);


parfor i = 1:length(pid_idx)
    % Define log file for the current patient
    log_filename = fullfile(plots_folder, ['patient_log_PC_' num2str(pid_idx(i)) '.txt']);
    fid = fopen(log_filename, 'w'); % Open the log file in write mode
    
    try
        index = find(Rid == Rid(pid_idx(i)));
        if length(index) > 2 
            fprintf(fid, '======index(end) of PID is %d\n', index(end));
            Num_pts = [Num_pts; length(index)];
            
            uv_Lp0 = zeros(136, 1);
            para = para0_all(i, :)'; % (137, 1)
            u0_para = u0_all(i, :)'; % (68, 1)
            w = w_all(i);
            patient_record = [];  % param for each patient for each random initial guess
            grad_vec = [];
            tspan = [50; Age(index(1:end-1)); 100];
            yy = Abeta(index(1:end-1), :)'; % (68, num_t)
            A_local = A_const.Value;
            A_local = double(A_local);
            Lp_local = Lp_const.Value;
            Lp_local = double(Lp_local);
            w_su = 1000; % weight for sparsity of u.
            w_sv = 1000; % weight for sparsity of v. 

            for k = 1:4
                fprintf(fid, 'Iteration %d: w_su=%f\n', k, w_su);
                options = optimoptions('fmincon', 'MaxIterations', 1e7);
                [uv_Lp, fval, exitflag, output, lambda_out, grad, ~] = fmincon(@(uv_Lp) loss_fun_normal_Abeta_68D_Lp_u(uv_Lp, para, u0_para, tspan, yy, double(Lp), w, w_u, w_su, w_v, w_sv), uv_Lp0, [], [], [], [], lb, ub, @(x) constraintFunction(x, A_local, 68), options);
                
                uv_Lp = double(uv_Lp); % (136, 1)
                grad_vec = grad(:);
                fprintf(fid, 'Exitflag: %d\n', exitflag);
                
                tspan1 = [50; Age(index(1:end)); 100];
                u_Lp = uv_Lp(1:68);
                v_Lp = uv_Lp(69:136);
                L = Lp + 0.5*[ones(1, 68) * v_Lp * diag(u_Lp) + ones(1, 68) * u_Lp * diag(v_Lp)] - 0.5 * (u_Lp * v_Lp' + v_Lp * u_Lp');
                fprintf(fid, 'Max eigenvalue of L: %f\n', max(eig(L)));
                fprintf(fid, 'Min eigenvalue of L: %f\n', min(eig(L)));
                
                [t, y] = ode45(@(t, y) ODEmodel_original_68D(para(1), para(2:69), para(70:137), double(L), y, t), tspan1, u0_para);
                err = y(ismember(t, Age(index(1:end))), :)' - Abeta(index(1:end), :)';
                pred_err_lsq = sqrt(sum(err.^2, 1)) ./ sqrt(sum(Abeta(index(1:end), :).^2, 2))'; % sum over rows, keep column
                
                fit_loss = sqrt(sum(pred_err_lsq(1:end-1).^2) / (length(index) - 1));
                err_T = sqrt(sum((y(ismember(t, 100), :)' - 1).^2));
                fprintf(fid, 'Fit loss: %f\n', fit_loss);
                fprintf(fid, 'Error at T=100: %f\n', err_T);
                
                patient_record = [patient_record, [index(end); w; w_su; w_sv; w_u; w_v; uv_Lp; fval; fit_loss; pred_err_lsq(end); exitflag; pred_err_lsq'; err_T]];
                w_su = w_su / 10;
                w_sv = w_sv / 10;
            end
            
            % Save patient records
            data_filename = fullfile(plots_folder, ['patient_record' num2str(index(end)) '.csv']);
            writematrix(patient_record, data_filename);

            data_filename = fullfile(plots_folder, ['patient_record_grad' num2str(index(end)) '.csv']);
            writematrix(grad_vec, data_filename);
        end
    catch ME
        % Log error information
        fprintf(fid, 'Error occurred: %s\n', ME.message);
    end
    
    % Close the log file
    fclose(fid);
end
