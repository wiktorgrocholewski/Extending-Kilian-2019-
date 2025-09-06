function [CI05,CI16,CI84,CI95] = fixedreg_boot(X,A_hat,u_hatf,u_hats,sigma_hat_uf,sigma_hat_us,p,Theta_hat,D0_hat,B,H,CI_method)
% Bootstrap Confidence Intervals for Structural IRFs with possible break
% Inputs:
% - X:        VAR regressors
% - A_hat:    Estimated VAR coefficients
% - u_hatf:   Residuals from first regime
% - u_hats:   Residuals from second regime
% - sigma_hat_uf / us: covariance estimates for each regime
% - p:        VAR lag length
% - Theta_hat: Point IRF estimate (K x K x H+1)
% - D0_hat:   Structural impact matrix
% - B:        Number of bootstrap replications
% - H:        Horizon for IRFs
% - CI_method: 'naive' or 't-percentile'

[T,Kp] = size(X);     
K = size(A_hat,1);    

A_boot = zeros(K,Kp,B);
Theta_boot = zeros(K,K,H+1,B);
t_stat = zeros(K,K,H+1,B);

for i = 1:B
    
    fprintf('Iteration %d out of %d\n', i, B);
    % Bootstrap residuals
    idxf = randi(size(u_hatf,1), size(u_hatf,1), 1);
    idxs = randi(size(u_hats,1), size(u_hats,1), 1);
    u_bootf = u_hatf(idxf,:) - mean(u_hatf(idxf,:),1);
    u_boots = u_hats(idxs,:) - mean(u_hats(idxs,:),1);
    u_boot = [u_bootf; u_boots];

    % Bootstrap data
    y_boot = (A_hat * X')' + u_boot;

    % VAR coefficients
    A_tmp = (X' * X) \ (X' * y_boot); 
    A_boot(:,:,i) = A_tmp';            
    Theta_b = irf(squeeze(A_boot(:,:,i)), D0_hat, H);
    Theta_boot(:,:,:,i) = Theta_b;

    if strcmp(CI_method, 't-percentile')
        % Residuals from boot VAR
        u_hatc = y_boot - X * A_boot(:,:,i)';
        u_hatcf = u_hatc(1:size(u_hatf,1), :);
        u_hatcs = u_hatc(size(u_hatf,1)+1:end, :);

        % Re-estimate impact matrix D0 for inner draw
        Sigma_hat_uc1 = cov(u_hatcf);
        Sigma_hat_uc2 = cov(u_hatcs);
        [D0_hatc, ~] = decompose_hetero(Sigma_hat_uc1, Sigma_hat_uc2); 

        % Inner bootstrap for std estimation
        Theta_bootc = zeros(K,K,H+1,B);
        for ii = 1:B
            idxfc = randi(size(u_hatcf,1), size(u_hatcf,1), 1);
            idxsc = randi(size(u_hatcs,1), size(u_hatcs,1), 1);
            u_bootcf = u_hatcf(idxfc,:) - mean(u_hatcf(idxfc,:),1);
            u_bootcs = u_hatcs(idxsc,:) - mean(u_hatcs(idxsc,:),1);
            u_bootc = [u_bootcf; u_bootcs];

            y_bootc = X * A_boot(:,:,i)' + u_bootc;
            A_tmpc = (X' * X) \ (X' * y_bootc);
            Theta_bootc(:,:,:,ii) = irf(A_tmpc', D0_hatc, H);
        end

        % Compute t-statistic (note: centered around Theta_hat)
        t_stat(:,:,:,i) = (Theta_boot(:,:,:,i) - Theta_hat) ./ std(Theta_bootc, 0, 4);
    end
end

% Compute confidence intervals
if strcmp(CI_method, 'naive')
    CI05 = squeeze(quantile(Theta_boot,0.05,4));
    CI16 = squeeze(quantile(Theta_boot,0.16,4));
    CI84 = squeeze(quantile(Theta_boot,0.84,4));
    CI95 = squeeze(quantile(Theta_boot,0.95,4));

elseif strcmp(CI_method, 't-percentile')
    t_stat(isnan(t_stat)) = 0;

    q05 = squeeze(quantile(t_stat, 0.05, 4));
    q16 = squeeze(quantile(t_stat, 0.16, 4));
    q84 = squeeze(quantile(t_stat, 0.84, 4));
    q95 = squeeze(quantile(t_stat, 0.95, 4));

    std_outer = std(Theta_boot, 0, 4);

    CI05 = Theta_hat + q05 .* std_outer;
    CI16 = Theta_hat + q16 .* std_outer;
    CI84 = Theta_hat + q84 .* std_outer;
    CI95 = Theta_hat + q95 .* std_outer;
else
    error('Specify CI_method correctly: "naive" or "t-percentile"');
end
end