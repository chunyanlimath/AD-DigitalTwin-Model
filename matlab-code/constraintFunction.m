function [c, ceq] = constraintFunction(x, A_p, n)
    % Extract u and v from x
    u = x(1:n);       % First n elements correspond to u
    v = x(n+1:end);   % Last n elements correspond to v
    
    % Compute A
    A = A_p + 0.5 * (u * v' + v * u') - diag(u .* v);
    
    % Inequality constraints
    c = [-A(:);  % Ensures 0 <= A_ij
         A(:) - 1]; % Ensures A_ij <= 1
    
    % No equality constraints
    ceq = [];
end


% function [c, ceq] = constraintFunction(x, A_p, n)
% % Extract u and v from x
%     u = x(1:n);       % First n elements correspond to u
%     v = x(n+1:end);   % Last n elements correspond to v
%     
%     % Compute A
%     A = A_p + 0.5 * (u * v' + v * u') - diag(u .* v);  % Calculate A matrix
%     
%     % Initialize constraints
%     c = zeros(2 * n * n, 1); % Preallocate for 2 constraints per pair (i, j)
%     ceq = [];  % No equality constraints
%     
%     % Populate constraints
%     idx = 1; % Index for constraint vector
%     for i = 1:n
%         for j = 1:n
%             % Constraints for each A_ij
%             c(idx) = -A(i, j);           % For 0 <= A_ij
%             c(idx + 1) = A(i, j) - 1;    % For A_ij <= 1
%             idx = idx + 2;
%         end
%     end
% end
%   