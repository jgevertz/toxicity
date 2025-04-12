%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Create Kaplan-Meier type curves of time to progression across VPs at  % 
% optimal protocol (by averaging across VPs) in the following cases:    %
% 1) Optimal adaptive protocol: Selected by looking across tumor size   %
%    thresholds and taking the on/off combination with highest average  %
%    TTP. (Data read in from Output/output_VPs.mat)                     %
% 2) Optimal dailty + toxicity protocol: Selected by looking across     %
%    toxicity thresholds and taking the on/off combination with highest %
%    average TTP. (Data read in from Output/output_VPs.mat)             %
% 3) Optimal adaptive + toxicity protocol: Selected by looking across   %
%    tumor size AND toxicity thresholds and taking the on/off           %
%    combination with highest average TTP. (Data read in from           %
%    Output/output_VPs_4D.mat)                                          %
% 4) Daily protocol: response across VPs. There are no thresholds to    %
%    optimize here.  (Data read in from Output/output_VPs.mat)          %
%                                                                       %
% Author: Jana Gevertz. Last update:2/24/2025                           %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all; 
num_VPs = 100; % Number of virtual patients 
num_pts = 13; % Discretization of protocol thresholds (use 11)
[q, ICs] = set_parameters_ICs();

%% Protocol settings
tspan = 0:0.1:100;
nrun = max(tspan);           % number of therapy cycles for continuous
baseline = (ICs(1) + ICs(2));
protocol_string = {'', 'Daily', 'Adaptive', 'Daily + Toxicity', 'Adaptive + Toxicity'};

%% KM curve settings and data
bin_length = 4; % Length of bin for plotting KM curve (in days)
KM_bins = 0:bin_length:nrun; 
load Output_VPs/output_VPs.mat % Load VP data

%% Daily only
daily_only = agg_tfail_tox_irrelevant(:,2);
%% Make KM plot for all VPs at protocol (Nhi_maxIdx, Nlo_max_idx)
KM_daily = zeros(num_VPs,1);
for i = 1:num_VPs
    KM_daily(i) = daily_only(i); 
    if KM_daily(i)==150
        KM_daily(i) = nan;
    end
end
KM_proportion_daily = size(KM_bins);
count_does_progress = 0;
for i = 2:length(KM_bins)
    %fprintf('From %d to %d:\n',KM_bins(i-1),KM_bins(i));
    A = find( (KM_daily<=KM_bins(i)) & (KM_daily>KM_bins(i-1)) );
    count_does_progress = count_does_progress + length(A);
    %fprintf('Found %d in bin #%d, making total = %d\n',length(A),i-1,count_does_progress)
    KM_proportion_daily(i-1) = 1 - (count_does_progress/num_VPs);
end
KM_proportion_daily(i) = KM_proportion_daily(i-1); % anyone who survived to 100
figure; 
stairs(KM_bins, KM_proportion_daily,'LineWidth', 2);
xlabel('Time (days)','FontSize',14);
ylabel('Proportion NOT Progressed','FontSize',14);
ylim([0,1])
title('VP Variability at Daily Protocol','FontSize',16);
% fname_fig = ['Output_VPs/KM_optimal_' protocol_string{2}];
% saveas(gcf,[fname_fig,'.fig'])
% saveas(gcf,[fname_fig,'.png'])


%% Adaptive only
adaptive_only = agg_tfail_sweep_adapt{3};
% adaptive_only(:,:,i) gives response in protocol space for VP #i
mean_adaptive_only = mean_tfail_sweep_adapt{3}; 
max_val_adapt = max(mean_adaptive_only,[],'all');
[Nhi_maxIdx,Nlo_maxIdx] = find(mean_adaptive_only==max_val_adapt);
fprintf('Adaptive only has max avg TTP of %f at relative Nhi = %f, relative Nlo = %f\n',...
    max_val_adapt,Nhi_vec(Nhi_maxIdx)/baseline,Nlo_vec(Nlo_maxIdx)/baseline);

%% Make KM plot for all VPs at protocol (Nhi_maxIdx, Nlo_max_idx)
KM_adaptive = zeros(num_VPs,1);
for i = 1:num_VPs
    KM_adaptive(i) = adaptive_only(Nhi_maxIdx,Nhi_maxIdx,i); 
    if KM_adaptive(i)==150
        KM_adaptive(i) = nan;
    end
end
KM_proportion_adaptive = size(KM_bins);
count_does_progress = 0;
for i = 2:length(KM_bins)
    A = find( (KM_adaptive<=KM_bins(i)) & (KM_adaptive>KM_bins(i-1)) );
    count_does_progress = count_does_progress + length(A);
    KM_proportion_adaptive(i-1) = 1 - (count_does_progress/num_VPs);
end
KM_proportion_adaptive(i) = KM_proportion_adaptive(i-1); % anyone who survived to 100
figure; 
stairs(KM_bins, KM_proportion_adaptive,'LineWidth', 2);
xlabel('Time (days)','FontSize',14);
ylabel('Proportion NOT Progressed','FontSize',14);
ylim([0,1])
title(['VP Variability at Optimal ' protocol_string{3} ' Protocol'],...
    'FontSize',16);
subtitle(['Max(avg TTP) = ' num2str(max_val_adapt) ' at relative Nhi = ' ...
    num2str(Nhi_vec(Nhi_maxIdx)/baseline) ', relative Nlo = ' ...
    num2str(Nlo_vec(Nlo_maxIdx)/baseline)],'FontSize',14);
% fname_fig = ['Output_VPs/KM_optimal_' protocol_string{3}];
% saveas(gcf,[fname_fig,'.fig'])
% saveas(gcf,[fname_fig,'.png'])


%% Toxicity only
toxicity_only = agg_tfail_sweep_tox{4};
% toxicity_only(:,:,i) gives response in protocol space for VP #i
mean_toxicity_only = mean_tfail_sweep_tox{4}; 
max_val_tox = max(mean_toxicity_only,[],'all');
[Tlo_maxIdx,Thi_maxIdx] = find(mean_toxicity_only==max_val_tox);
fprintf('Dailty+tox has max avg TTP of %f at Tlo = %f, Thi = %f\n',...
    max_val_tox,Tlo_vec(Tlo_maxIdx),Thi_vec(Thi_maxIdx));

%% Make KM plot for all VPs at protocol (Tlo_maxIdx, Thi_max_idx)
KM_toxicity = zeros(num_VPs,1);
for i = 1:num_VPs
    KM_toxicity(i) = toxicity_only(Tlo_maxIdx,Thi_maxIdx,i); 
    if KM_toxicity(i)==150
        KM_toxicity(i) = nan;
    end
end
KM_proportion_toxicity = size(KM_bins);
count_does_progress = 0;
for i = 2:length(KM_bins)
    A = find( (KM_toxicity<=KM_bins(i)) & (KM_toxicity>KM_bins(i-1)) );
    count_does_progress = count_does_progress + length(A);
    KM_proportion_toxicity(i-1) = 1 - (count_does_progress/num_VPs);
end
KM_proportion_toxicity(i) = KM_proportion_toxicity(i-1); % anyone who survived to 100
figure; 
stairs(KM_bins, KM_proportion_toxicity,'LineWidth', 2);
xlabel('Time (days)','FontSize',14);
ylabel('Proportion NOT Progressed','FontSize',14);
ylim([0,1])
title(['VP Variability at Optimal ' protocol_string{4} ' Protocol'],...
    'FontSize',16);
subtitle(['Max(avg TTP) = ' num2str(max_val_tox) ' at Tlo = ' ...
    num2str(Tlo_vec(Tlo_maxIdx)) ', Thi = ' num2str(Thi_vec(Thi_maxIdx))],...
    'FontSize',14);
% fname_fig = ['Output_VPs/KM_optimal_' protocol_string{4}];
% saveas(gcf,[fname_fig,'.fig'])
% saveas(gcf,[fname_fig,'.png'])

%% Now the 4D protocol sweep
load Output_VPs/output_VPs_4D.mat % Load VP data
fprintf('Adaptive+tox has max avg TTP of %f with:\n',max_val_all);
fprintf('\tToxicity thresholds: Tlo = %f, Thi = %f\n',...
    Tlo_vec(Tlo_maxIdx),Thi_vec(Thi_maxIdx));
fprintf('\tAdaptive thresholds: Relative Nhi = %f, relative Nlo = %f\n',...
    Nhi_vec(Nhi_maxIdx)/baseline,Nlo_vec(Nlo_maxIdx)/baseline);

%% Make KM plot for all VPs at protocol (Tlo_maxIdx, Thi_max_idx)
KM_adapt_tox = zeros(num_VPs,1);
for i = 1:num_VPs
    KM_adapt_tox(i) = agg_tfail_sweep_all(Tlo_maxIdx,Thi_maxIdx,...
        Nhi_maxIdx,Nlo_maxIdx,i); 
    if KM_adapt_tox(i)==150
        KM_adapt_tox(i) = nan;
    end
end
KM_proportion_adapt_tox = size(KM_bins);
count_does_progress = 0;
for i = 2:length(KM_bins)
    A = find( (KM_adapt_tox<=KM_bins(i)) & (KM_adapt_tox>KM_bins(i-1)) );
    count_does_progress = count_does_progress + length(A);
    KM_proportion_adapt_tox(i-1) = 1 - (count_does_progress/num_VPs);
end
KM_proportion_adapt_tox(i) = KM_proportion_adapt_tox(i-1); % anyone who survived to 100
figure; 
stairs(KM_bins, KM_proportion_adapt_tox,'LineWidth', 2);
xlabel('Time (days)','FontSize',14);
ylabel('Proportion NOT Progressed','FontSize',14);
ylim([0,1])
title(['VP Variability at Optimal ' protocol_string{5} ' Protocol'],...
    'FontSize',16);
subtitle(['Max(avg TTP) = ' num2str(max_val_tox) ' at Tlo = ' ...
    num2str(Tlo_vec(Tlo_maxIdx)) ', Thi = ' num2str(Thi_vec(Thi_maxIdx)) ...
    ', rel. Nhi = ' num2str(Nhi_vec(Nhi_maxIdx)/baseline) ', rel. Nlo = ' ...
    num2str(Nlo_vec(Nlo_maxIdx)/baseline)],'FontSize',11);
% fname_fig = ['Output_VPs/KM_optimal_' protocol_string{5}];
% saveas(gcf,[fname_fig,'.fig'])
% saveas(gcf,[fname_fig,'.png'])


%% All together now
figure; 
stairs(KM_bins, KM_proportion_daily,'LineWidth', 2);
hold on;
stairs(KM_bins, KM_proportion_adaptive,'LineWidth', 2);
stairs(KM_bins, KM_proportion_toxicity,'LineWidth', 2);
stairs(KM_bins, KM_proportion_adapt_tox,'LineWidth',2)
hold off;
xlabel('Time (days)','FontSize',14);
ylabel('Proportion NOT Progressed','FontSize',14);
ylim([0,1])
title('VP Variability at Optimal Protocol','FontSize',16);
legend(protocol_string{2},protocol_string{3},protocol_string{4},...
    protocol_string{5},'Location','bestoutside','FontSize',14); 
legend boxoff
fname_fig = 'Output_VPs/KM_optimal_all_withDaily';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])


%% Updated set_parameters_ICs function
function [q, ICs] = set_parameters_ICs()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%     Tumor Growth/Death Params     %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    q.alphaS = 1;           % per unit time, sensitive cell growth rate
    q.epsilon = 0.4;          % Default epsilon for initial setup
    q.alphaR = q.epsilon * q.alphaS; % Resistant cell growth rate defined via epsilon

    q.K = 100;              % Carrying capacity
    q.beta = 2.4;           % Competitive advantage/disadvantage of sensitive  
                            % cells over resistant cells
                            % beta > 1 => advantage
                            % beta < 1 => disadvantage
                            % beta = 1 => neutral
    q.delta = 1;            % per unit time per drug, sensitive cell death rate

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%     Drug Params     %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    q.lambda = log(2);      % Assuming drug half life is 1 unit time
    q.mu = 1;               % time scale of tox onset
    q.gamma = 0.4;          % Rate at which drug Tox declines

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  Initial Condtions  %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S0 = 49.5;     % percent of capacity that is sensitive cells
    R0 = 0.5;      % percent of capacity that is resistant cells
    D0 = 0;        % Initial drug concentration
    T0 = 0;        % Initial toxicity level

    ICs = [S0, R0, D0, T0];
end

