function plot_irf(Theta_hat, CI05, CI16, CI84, CI95, shock)
% Plots the estimated IRF with 95% and 68% Confidence Intervals
% and a checkered grid
% INPUT
% Theta_hat : K x K x H+1 Estimated Structural IRF
% CI05/CI95 : 90% CI bounds
% CI16/CI84 : 68% CI bounds
% shock     : Index of shock to plot responses to

K = size(Theta_hat,1);
H = size(Theta_hat,3)-1;
h = 0:H;

varNames_paper = {'Real oil price','World oil production','World oil inventories','Global real activity','U.S. industrial production','U.S. CPI'};

% Colors
blue90 = [0.7, 0.9, 1.0];      % Light blue for 90% CI
blue68 = [0.4, 0.65, 0.85];    % Darker blue for 68% CI
alpha90 = 0.4;
alpha68 = 0.6;
lineColor = [0.1, 0.1, 0.1]; % Near black

figure()
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

for j = 1:K
    subplot(2, K/2, j)
    hold on

    % 90% CI patch
    xQ1 = squeeze(CI05(j,shock,:));
    xQ2 = squeeze(CI95(j,shock,:));
    xx = [h, fliplr(h)];
    yy = [xQ1' fliplr(xQ2')];
    patch(xx, yy, blue90, 'EdgeColor', 'none', 'FaceAlpha', 1);

    % 68% CI patch
    yQ1 = squeeze(CI16(j,shock,:));
    yQ2 = squeeze(CI84(j,shock,:));
    yy = [yQ1' fliplr(yQ2')];
    patch(xx, yy, blue68, 'EdgeColor', 'none', 'FaceAlpha', 1);

    % Zero line
    plot(h, zeros(size(h)), 'k-', 'LineWidth', 1);

    % IRF line
    irf = squeeze(Theta_hat(j,shock,:));
    plot(h, irf, '-', 'Color', lineColor, 'LineWidth', 2.5);

    % Grid styling (checkered background)
    grid on
    ax = gca;
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';
    ax.GridAlpha = 0.2;
    ax.GridColor = [0.5, 0.5, 0.5];

    % Formatting
    set(gca, 'fontsize', 20);
    xlim([0 H]);
    xlabel('Months');
    ylabel('%');
    title(varNames_paper{j});

    if j == 1
        legend({'90% CI', '68% CI', 'Zero', 'IRF'}, 'Location', 'Best');
    end

    hold off
end

end