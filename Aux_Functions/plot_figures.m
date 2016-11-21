% Script for Plotting Data
%
% See COPYRIGHT for license information 
%

clear; close all; clc

addpath SMC_Code
addpath Model
addpath Aux_Functions
addpath Results

% Load Parameters
params = struct;
params.alw = 1;      % AxesLineWidth
params.fsz = 16;     % Fontsize
params.lw  = 2;      % LineWidth
params.msz = 4;      % MarkerSize

input_text   = {'$H$','$P_{o}$', '$F_s$','$\theta$','$L_{sp}$','$q$','$RD$','$L_{a}$', '$C_{dps}$'};
output_text  = {'$\tau_{tot}$','$P_{tot}$','$A_{sa}$'};
input_label  = {'H','Po', 'Fs','theta','Lsp','q','RD','La', 'Cdps'};
output_label = {'Tto','Ptot','Asa'};

mkdir('Figures');

% Determine the number of inputs and outputs
n_inputs  = length(input_text);
n_outputs = length(output_text);

% Load data file
load(['Results/post_SMC_data'])

ax_lims = {[0, 0.01, 0.02, 0.03, 0.04], [900, 1000, 1100, 1200, 1300, 1400], [10, 11, 12, 13, 14]};

% Plot data for original model using KDE
ref_axes = matrix_plot(ref_model_mc_full_output, output_text, params);
print('-depsc',['Figures/model_dist_kde_ref']);

for i=1:length(lambda_vect)

	% Plot data using Histogram
	mc_gauss_comp(model_mc_full_output{i}, opt_model_gaussian{i,1}, opt_model_gaussian{i,2}, ...
				  ps.ref_mean, ps.ref_cov, output_text, params);
	print('-dpng',['Figures/model_dist_histogram_' lambda_str{i}])

	% Plot data using KDE
	matrix_plot_update(model_mc_full_output{i}, output_text, params, ref_axes);
	print('-depsc',['Figures/model_dist_kde_' lambda_str{i}]);

	n_outputs = size(model_mc_full_output{i},2);
	for j=1:n_outputs
		figure;
		hold on
		grid on
		scatter(ref_model_mc_full_output(:,j), model_mc_full_output{i}(:,j),'o','MarkerEdgeColor',[0, 0, 0.8],'MarkerEdgeAlpha',0.3)
		%scatter(ref_model_mc_full_output(:,j), model_mc_full_output{i}(:,j),'o','MarkerFaceColor', ...
		%		[0, 0, 0.8],'MarkerEdgeColor',[0, 0, 0.8],'MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.3)
		plot([ax_lims{j}(1), ax_lims{j}(end)],[ax_lims{j}(1), ax_lims{j}(end)],'-r','LineWidth', params.lw);
		xlabel([output_text{j} ' - Reference Model'],'interpreter','latex','FontSize',params.fsz)
		ylabel([output_text{j} ' - Decoupled Model'],'interpreter','latex','FontSize',params.fsz)
		set(gca,'LabelFontSizeMultiplier',1.2,'FontSize', params.fsz, 'LineWidth', params.alw);
		set(gca,'Ytick',ax_lims{j},'YLim',[ax_lims{j}(1), ax_lims{j}(end)])
		set(gca,'Xtick',ax_lims{j},'XLim',[ax_lims{j}(1), ax_lims{j}(end)])
		hold off
		iptsetpref('ImshowBorder','tight');
		print('-dpng',['Figures/qoi_comp_model_' lambda_str{i} '_' output_label{j}])
		print('-depsc',['Figures/qoi_comp_model_' lambda_str{i} '_' output_label{j}])
	end

end

figure;
semilogx(lambda_vect, opt_models_kl_div, '-o', 'linewidth',2);
hold on
grid on
xlabel('$\lambda$','interpreter','latex','FontSize',14)
ylabel('KL Divergence','interpreter','latex','FontSize',14)
hold off
print('-dpng','Figures/kl_divergence_vs_lambda')

figure;
plot(opt_models_coupls, opt_models_kl_div, '-o', 'linewidth',2);
hold on
grid on
xlabel('Retained Connections','interpreter','latex','FontSize',14)
ylabel('KL Divergence','interpreter','latex','FontSize',14)
hold off
print('-dpng','Figures/kl_divergence_vs_connections')

% -- END OF FILE --