
data = readmatrix('With_age_UCBERKELEYAV45_04_26_22.csv');
Rid = data(:, 1);  % Research participant identifier, unique for every patient

final_1d=readmatrix('plots_Abeta_scalar/fmincon_nolsq_w_0_00015625Lp/summary_best_w_patient_record.xlsx');
pid_idx = final_1d(:, 1);
record_all = zeros(471, 153);
 
plots_folder1 = 'plots_Abeta_68D_u0/fmincon_nolsq_w_Lp_68D_regL1_u0/';   
data_filename = fullfile(plots_folder1, 'summary_patient_record.xlsx');
results_all = readmatrix(data_filename);
fit_68D_u0 = results_all(:, 4);
pred_68D_u0 = results_all(:, 5);
plots_folder = 'plots_Abeta_Lp_uu/fmincon_nolsq_68D_u0_Lp_lowrank_u_wu1000_ws1000-cluster/';
para_u_Lp_all = [];
fit_Lp_all = zeros(471, 1);
pred_Lp_all = zeros(471, 1);
 

for i = 1:length(pid_idx)
    index = find(Rid == Rid(pid_idx(i)));
    disp(['----- pid idx=', num2str(index(end))]);
    

    % Read patient record from the two folders
    data_filename = fullfile(plots_folder, ['patient_record' num2str(index(end)) '.csv']);
    record = readmatrix(data_filename); % (153, 6)

    % Extract fit_loss and pred values
    fit_loss = record(144, :);
    pred = record(145, :)';

    % Find the best fit and corresponding index
    [best_fit, best_column_idx] = min(fit_loss);

    % Check if pred of the best_column_idx is not smaller than pred_68D_u0(i)
    if pred(best_column_idx) - pred_68D_u0(i) > 0
        best_column_idx = nan;  % Default to w_s = 0 if condition is not met
        record_all(i, 1:147 + length(index)) = [index(end), nan(1, 146 + length(index))];
        fit_Lp_all(i) = fit_68D_u0(i);
        pred_Lp_all(i) = pred_68D_u0(i);
        para_u_Lp_all(i, :) = [Rid(index(end)), zeros(136, 1)'];
        disp(['Go back to alpha=0,' num2str(i)]);
        % record_all(i, 1:length(record_u0)) = record_u0;
    else
        disp(['Best column idx=' num2str(best_column_idx)]);
        record_all(i, 1:147 + length(index)) = record(:, best_column_idx)';
        fit_Lp_all(i) = fit_loss(best_column_idx);
        pred_Lp_all(i) = pred(best_column_idx);
        para_u_Lp_all(i, :) = [Rid(index(end)), record(7:142, best_column_idx)'];
    end
end
% Store results in record_all
Lp_fit_pred = [fit_Lp_all, pred_Lp_all];
data_filename = fullfile(plots_folder, 'collect_best_w_s_fit_all_Lp_fit_pred_loss.xlsx');
writematrix(Lp_fit_pred, data_filename);

% Save the result to a file
data_filename = fullfile(plots_folder, 'collect_best_Lp_uu.xlsx');
writematrix(para_u_Lp_all, data_filename);

% Save the result to a file
data_filename = fullfile(plots_folder, 'collect_best_w_s_fit_all_Lp.xlsx');
writematrix(record_all, data_filename);

fit_loss=zeros(471, 1);
pred=zeros(471, 1);
for i=1:length(pid_idx)
    % index=find(Rid==Rid(pid_idx(i)));
    fit_loss(i)=record_all(i, 144);
    pred(i)=record_all(i, 145);
end
all=[fit_loss, pred];
find(fit_Lp_all-fit_68D_u0>1e-4)
find(pred-pred_68D_u0>1e-4)
