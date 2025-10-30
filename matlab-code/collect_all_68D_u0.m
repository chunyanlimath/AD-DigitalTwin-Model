
data = readmatrix('With_age_UCBERKELEYAV45_04_26_22.csv');
Rid = data(:, 1);
PID = unique(Rid);

plots_folder = 'plots_Abeta_68D_u0/fmincon_nolsq_w_Lp_68D_regL1_u0/';  %_homotopy_w

results_all=zeros(471, 80); 
results_all_u0 = zeros(471, 68);
k = 0;

final_1d=readmatrix('scalar_baseline_result.xlsx');
pid_idx=final_1d(:, 1);

for i = 1:length(pid_idx)
    index=find(Rid==Rid(pid_idx(i)));
    if length(index)>2
        k = k+1;
        fit_all = [];
        %  patient_record = [patient_record; [index(end); fval; fit_loss; pred_err_lsq(end); exitflag; pred_err_lsq'; err_T; u0_para]'];
               
        data_filename = fullfile(plots_folder, ['patient_record' num2str(index(end)) '.csv']);
        record = readmatrix(data_filename);
        results_all(k, 1:length(record))=record;
        %patient_record = [patient_record; [index(end); fval; fit_loss; pred_err_lsq(end); exitflag; pred_err_lsq'; err_T; u0_para]'];
        data_filename = fullfile(plots_folder, ['patient_record_u0' num2str(index(end)) '.csv']);
        record_u0 = readmatrix(data_filename); % a column vector
        results_all_u0(k, :)=record_u0'; 
    end
end


data_filename = fullfile(plots_folder, 'summary_patient_record.xlsx');
writematrix(results_all, data_filename);

data_filename = fullfile(plots_folder, 'summary_patient_record_u0.xlsx');
writematrix(results_all_u0, data_filename);

