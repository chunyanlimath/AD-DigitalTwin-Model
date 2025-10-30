% load scalar baseline model results 

% result scalar baseline model 
% header=[PID	w	D	lambda	K u0 fit_loss pred_loss	fval 	grad_D	grad_lambda	grad_K grad_u0	exitflag	t1	t2	t3	t=100];
% final_1d = readmatrix('scalar_baseline_result.xlsx');
final_1d=readmatrix('plots_Abeta_scalar/fmincon_nolsq_w_0_00015625Lp/summary_best_w_patient_record.xlsx');
pid=final_1d(:, 1);
w=final_1d(:, 2);
D=final_1d(:, 3);
lambda=final_1d(:, 4);
K=final_1d(:, 5);
Abeta0=final_1d(:, 6);
fit_scalar=final_1d(:, 7);
pred_scalar=final_1d(:, 8);
fval=final_1d(:, 9);

load('Lp68.mat');

% ------- Load the required data --------------------------------------------------------------------
data = readmatrix('With_age_UCBERKELEYAV45_04_26_22.csv');
Rid = data(:, 1);
PID = unique(Rid);
Age = data(:, 4);

Abeta_ref = data(:, 12);
Abeta = data(:, 42:2:181)./ Abeta_ref;
% remove two unknown columns 
columns_with_all_nan = all(isnan(Abeta), 1); % Calculate the indices of columns that are entirely NaN
Abeta = Abeta(:, ~columns_with_all_nan); % Remove those columns
disp('Size of Abeta after removing columns with all NaN values:');
disp(size(Abeta)); % (3086, 68)

% normalization 
Abeta_Max = max(max(Abeta));
Abeta = Abeta / Abeta_Max;

% load the results with L1 (\epslion)
plots_folder = 'plots_Abeta_68D/fmincon_Lp_68D_regL1_vary_wreg';
data_filename = fullfile(plots_folder, 'collect_best_w_reg_fitloss_Sshape.xlsx');
result_68d = readmatrix(data_filename);
pred_68d=result_68d(:, 5);
fit_68d=result_68d(:, 4);
 
plots_folder = 'plots_Abeta_68D_u0/fmincon_nolsq_w_Lp_68D_regL1_u0';
fval_eps=zeros(471, 1);
pred_eps=zeros(471, 1);
fit_loss_eps=zeros(471, 1);
epslion_eps=[];
record_eps = zeros(471, 6+6+1);  % record collection of best w_reg result 
record_eps_para = zeros(471, 68);  % collection of best para for best w_reg

Num_pts = [];
for i=1:length(pid)
    index = find(Rid==Rid(pid(i)));
    Num_pts = [Num_pts; length(index)];
    if length(index)>2                                        
        data_filename = fullfile(plots_folder, ['patient_record' num2str(pid(i)) '_w0_00390625.csv']);
        [record_1]=xlsread(data_filename);
        % patient_record = [patient_record; [index(end); w1; fval; fit_loss; pred_err_lsq(end); exitflag; pred_err_lsq'; err_T; u0_para+u0]'];
                    
        fit_loss=record_1(1:11, 4);       
        pred_loss=record_1(1:11, 5);
        % select best w_reg based on fit_loss and pred_loss compared with
        % scalar baseline model
        [best_fit_eps, best_row_idx] = min(fit_loss);
        len=6+1+length(index);
        PIDmy=record_1(:, 1);
        index = find(Rid == Rid(PIDmy(1)));
            
        disp(['====== PID is ' num2str(index(end))]);            
        if pred_loss(best_row_idx, end) - pred_68d(i) < 1e-4
           
           record_eps(i, 1:len)=record_1(best_row_idx, 1:len); % if pred is improved, then choose this one!
           record_eps_para(i,:)=record_1(best_row_idx, len+1:end);
        else 
           best_row_idx=1;
           record_eps(i, 1:len)=record_1(1, 1:len); % if pred is not improved, then, choose baseline which set w_reg=1000
           record_eps_para(i,:)=record_1(best_row_idx, len+1:end);
        end
        pred_eps(i)=pred_loss(best_row_idx); 
        fit_loss_eps(i)=fit_loss(best_row_idx);
    end
end

data_filename = fullfile(plots_folder, 'summary_patient_record.xlsx');
writematrix(record_eps, data_filename);

data_filename = fullfile(plots_folder, 'summary_patient_record_u0.xlsx');
writematrix(record_eps_para, data_filename);

pred_err = pred_eps-pred_scalar;
sum(pred_err)
figure;
plot(pred_err, 'bo'); % Plot pred_err with blue circles
hold on;
plot(fit_loss_eps-fit_scalar, 'r+'); % Plot fit_loss_eps with red pluses
hold on;
yline(0, 'k-');
hold off;
legend('pred\_eps - pred\_scalar', 'fit\_eps - fit\_scalar'); % Add a legend to label the plots
title('Comparison of Prediction Error and Fit Loss with scalar'); % Optional: Add a title
data_filename = fullfile(plots_folder,'comparison_based_on_fitloss_u0_scalar.png');
saveas(gcf, data_filename); % Save the figure as a PNG file
 
pred_err = pred_eps-pred_68d;
sum(pred_err)
figure;
plot(pred_err, 'bo'); % Plot pred_err with blue circles
hold on;
plot(fit_loss_eps-fit_68d, 'r+'); % Plot fit_loss_eps with red pluses
hold on;
yline(0, 'k-');
hold off;
legend('pred\_eps - pred\_scalar', 'fit\_eps - fit\_scalar'); % Add a legend to label the plots
title('Comparison of Prediction Error and Fit Loss with scalar'); % Optional: Add a title
data_filename = fullfile(plots_folder,'comparison_based_on_fitloss_u0_68d.png');
saveas(gcf, data_filename); % Save the figure as a PNG file
 