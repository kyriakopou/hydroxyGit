figure();
x=linspace(0, 1);
mu = 0.7;
sigma = 0.2;
plot(x, normpdf(x, mu, sigma), 'LineWidth', 4);
hold on;
yLim = get(gca, 'ylim');
plot([mu, mu], yLim, 'LineStyle', '--',  'LineWidth', 4);
text(mu+0.05, 1, '$\mathbf{x}_{i-1}$', 'Interpreter', 'latex', 'FontSize', 25);

%remove axes ticks and set the essential ones
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 25);
set(gca,'ytick',{});
set(gca,'xtick',[0 1]);

%plot the sample
sample = 0.5;
plot(sample, normpdf(sample, mu, sigma), 'o', 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor',  [0.466   0.674   0.188])
text(sample-0.15, 1.3, '$\mathbf{y}_{i-1}$', 'Interpreter', 'latex', 'FontSize', 25);

print('~/Desktop/normal_prior_inter', '-dpdf');

%bounds for the gradients
a = -sample/6;
b = (1-sample)/6;
x = linspace(a, b);
mu_grad = 0.04;
sigma_grad = 0.02;

figure();
plot(x, normpdf(x, mu_grad, sigma_grad), 'LineWidth', 4);
yLim = get(gca, 'ylim');
xlim([x(1) x(end)]);
hold on;
plot([mu_grad, mu_grad], yLim, 'LineStyle', '--',  'LineWidth', 4);
text(mu_grad+0.01, 10, '$\mathbf{x}_{i}$', 'Interpreter', 'latex', 'FontSize', 25);


%remove axes ticks and set the essential ones
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 25);
set(gca,'ytick',{});
set(gca,'xtick', [x(1) x(end)]);
% set(gca,'XTickLabel', {'$a_i=\phi_a(\mathbf{v}_{i-1})$', '$b_i = \phi_b(\mathbf{v}_{i-1})$'}, 'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'$a_i$', '$b_i$'}, 'TickLabelInterpreter', 'latex');


%plot the sample
sample = 0.053;
plot(sample, normpdf(sample, mu_grad, sigma_grad), 'o', 'LineWidth', 2, 'MarkerSize', 15, 'MarkerFaceColor',  [0.466   0.674   0.188])
text(sample+0.01, normpdf(sample, mu_grad, sigma_grad)+1, '$\mathbf{y}_{i}$', 'Interpreter', 'latex', 'FontSize', 25);

print('~/Desktop/normal_prior_grad', '-dpdf');
close all;
