% in this code, we optimize the model w.r.t para=[D, lamda, K, Abeta0] are collection of 3 scalars and intial. 
% Now we vary w (homotopy of w) to searching for best w to reach the optimial value of objective function 
% value of w is optimized for each person. 
% in this code, we try to vary w start from a large value, say 0.01, reduce to a smaller value, say, 1e-4.
% L is L_p which is the mean of all Laplacian 

% load Lp matrix with size 68 * 68. 
load('Lp68.mat');
% ---------- load data ------------------ 
data = readmatrix('With_age_UCBERKELEYAV45_04_26_22.csv');
Rid=data(:,1); % research participant identifier which is unique for every patient
PID=unique(Rid); % one patient might have multiple fMRI scan or PET scan 
Age=data(:,4);

Abeta_ref = data(:, 12); % COMPOSITE_REF_SUVR 
Abeta=data(:,42:2:181)./Abeta_ref;  % Abeta info in Cortext region
columns_with_all_nan = all(isnan(Abeta), 1); % Calculate the indices of columns that are entirely NaN
Abeta = Abeta(:, ~columns_with_all_nan); % Remove those columns
disp('Size of Abeta after removing columns with all NaN values:');
disp(size(Abeta)); % (3086, 68)
% --------- normalization ----------------
Abeta_Max = max(max(Abeta))
Abeta = Abeta / Abeta_Max;

Accuracy=[];
best_fit=[];
Num_pts=[];  
exitflags = [];
rng(567093); % Set the seed for reproducibility
N_initial = 50;
D = rand(N_initial, 1); % D in [0, 1]
lambda =rand(N_initial, 1); % lambda in [0, 1]
K = 1.0 + rand(N_initial, 1);  % K in [1, 2]
 
plots_folder =  'plots_Abeta_scalar/fmincon_nolsq_w_0_00015625Lp';

% Create the folder if it doesn't exist
if ~exist(plots_folder, 'dir')
    mkdir(plots_folder);
end

parfor i=1:length(PID)
    index=find(Rid==PID(i));
    % Define log file for the current patient
    log_filename = fullfile(plots_folder, ['patient_log_' num2str(index(end)) '.txt']);
    fid = fopen(log_filename, 'w'); % Open the log file in write mode
    
    try
        if length(index)>2 
            disp(['======index(end) of PID is ' num2str(index(end))]);
            L = Lp;
            Num_pts=[Num_pts;length(index)];
            lb = [0.0; 0.0; 1.0; 0.0]; % (4, 1) user defined lower bounds for [D, lambda, K, Abeta0]
            ub = [1.0; 1.0; 2.0; 1.0]; % (4, 1) user defined
            patient_record_all = [];  % param for each patient for all random initial guess   

            for idx = 1:N_initial
                disp(['====== rand i.c. idx is ' num2str(idx)]);
                para0 = [D(idx); lambda(idx); K(idx); rand(1)]; % (4, 1)
                fprintf('the para0=[D, lambda, K, Abeta0] is %.4f, %.4f, %.4f, %.4f \n', para0(1), para0(2), para0(3), para0(4));
                patient_record = [];  % param for each patient for each random initial guess
                w=0.0025;
                for k = 1:4 % w=0.0025/2^{1, ..., 3}
                    disp(['w=', num2str(w)]);
                    tspan = [50; Age(index(1:end-1)); 100]; % (68, num_t-1 + 1)
                    yy = Abeta(index(1:end-1), :)'; % (68, num_t-1)
                    [para, fval, exitflag, output, lambda_out, grad, ~]=fmincon(@(para) loss_fun_normal_Abeta_scalar(para, tspan, yy, double(L), w),para0,[],[],[],[],lb, ub,[],optimset('MaxFunEvals',1e7));
                    para = double(para);
                   
                    grad_vec = grad(:);

                    disp(['exitflag=', num2str(exitflag)]);
                    tspan1 = [50; Age(index(1:end)); 100];
                    [t, y]=ode45(@(t,y) ODEmodel_original(para(1), para(2), para(3), double(L), y, t),tspan1, para(4)*ones(68, 1));
                    % y size = (num_t, 68)

                    err=y(ismember(t, Age(index(1:end))), :)' - Abeta(index(1:end), :)'; % (68, num_t)
                    pred_err_lsq = sqrt(sum(err.^2, 1)) ./ sqrt(sum(Abeta(index(1:end), :).^2, 2))'; % (1, num_t)
                    
                    fit_loss = sqrt(sum(pred_err_lsq(1:end-1).^2)/(length(index)-1)); 
                    err_T = sqrt(sum((y(ismember(t, 100), :)' - 1).^2));
                    patient_record=[patient_record, [index(end); w; para; fit_loss; pred_err_lsq(end); fval; grad_vec; exitflag; pred_err_lsq'; err_T]];
                    disp(['fmincon err and err for T=100 is ', num2str(pred_err_lsq) ' ' num2str(err_T)]);
                    disp(['fmincon para=[D, lambda, K, y0] is ', num2str(para')]);

                    % Reduce w by a factor of 1/2 and update para0
                    w = w/2; 
                    para0 = para; % homotopy idea for w
                end 
                patient_record_all = [patient_record_all; patient_record];
            end
            % Generate the folder and filename for saving the results
            w_str = sprintf('%.8f', w);   
            w_str = strrep(w_str, '.', '_');
 
            data_filename = fullfile(plots_folder, ['patient_record_grad' num2str(index(end)) '_rand_idx' num2str(idx) '_' w_str '.csv']); %homotopy_w
            writematrix(patient_record_all, data_filename);
            
        end
    catch ME
        % Log error information
        fprintf(fid, 'Error occurred: %s\n', ME.message);
    end
    
    % Close the log file
    fclose(fid);
end


