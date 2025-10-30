load('Lp68.mat');
% check the results of scalar baseline model 
% make sure that the value w is selected based on the fitting error instead of fval or prediction error
% ------- Load the required data -------------------------------------------------------------------
[data, txt] = xlsread('With_age_UCBERKELEYAV45_04_26_22.csv');
Rid = data(:, 1);
PID = unique(Rid);
Age = data(:, 4);

Abeta_ref = data(:, 12);
Abeta = data(:, 42:2:181)./ Abeta_ref;
columns_with_all_nan = all(isnan(Abeta), 1); % Calculate the indices of columns that are entirely NaN
Abeta = Abeta(:, ~columns_with_all_nan); % Remove those columns
disp('Size of Abeta after removing columns with all NaN values:');
disp(size(Abeta)); % (3086, 68)
%  normalization 
Abeta_Max = max(max(Abeta));
Abeta = Abeta/Abeta_Max;

best_fit = [];
Num_pts = [];

plots_folder = 'plots_Abeta_scalar/fmincon_nolsq_w_0_00015625Lp/';
data_filename = fullfile(plots_folder, 'summary_patient_record_original.xlsx');
[record] = xlsread(data_filename);
% [PID, w, D, lambda, K, u0, fit_loss, pred_loss, fval, grad_D, grad_lambda, grad_K, grad_u0, exitflag, t1, t2, ...]
record_all = zeros(471, 21);
pred_scalar = zeros(471, 1);
% [index(end); w; para; fit_loss; pred_loss; fval; grad_vec; exitflag; pred_err_lsq'; err_T]
for i = 1:length(PID)
    index=find(Rid==PID(i));
    if length(index)>2
        idx_start = (i-1) * 4 + 1;
        idx_end = i * 4;
        pid = record(idx_start:idx_end, 1);
        disp(['pid= ' num2str(mean(pid))]);
        fit_loss = record(idx_start:idx_end, 7);
        pred_scalar = record(idx_start:idx_end, 8); 
        fit_loss_copy=fit_loss;
        attempts = 0;
        foundNonDecreasing = false;
        index = find(Rid == Rid(pid(1)));
        disp(['====== PID is ' num2str(index(end))]);
        L = Lp;

        tspan1 = [50; Age(index(1:end)); 100];
                 
        % best_row_idx=nan;% initial 
        while attempts < 4 && ~foundNonDecreasing
            [best_fit_eps, best_row_idx] = min(fit_loss);
            para = record(idx_start+best_row_idx-1, 3:6);
            [t, y] = ode45(@(t, y) ODEmodel_original(para(1), para(2), para(3), double(L), y, t), tspan1 , para(4)*ones(68, 1));
            isNonDecreasing = all(diff(y) >= 0); % logical array
            sum(isNonDecreasing)
            if sum(isNonDecreasing) > 65
                record_all(i, :) = record(idx_start + best_row_idx - 1, :);
                foundNonDecreasing = true;
            else
                fit_loss(best_row_idx, :) = nan;
                attempts = attempts + 1;
            end
        end

        if ~foundNonDecreasing
             disp(['No non-decreasing solution found for PID ' num2str(PID(i)) '. Fallback to the first row.']);
             record_all(i, :) = record(idx_start, :);
        end
    end
end

% Create a logical index for rows that are not all zeros
nonZeroRows = any(record_all ~= 0, 2);
% Remove rows that are all zeros
results_all_filtered = record_all(nonZeroRows, :);

data_filename = fullfile(plots_folder, ['summary_best_w_patient_record.xlsx']);
writematrix(results_all_filtered, data_filename);