%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Quantifying the impact of most sensitive parameters and treatment-    %
% related thresholds at the nominal model parameterization.             %
% Updated: 2/24/2025                                                    %
%                                                                       %
% Parameters swept over:                                                %
% - Competition parameter, beta                                         %
% - Resistant cell growth rate, alpha_R = epsilon*alpha_S, 0<eps<1      %
% For each parameterization, code determines the protocol (daily,       %
% adaptive, daily+tox, adaptive+tox) that maximizes time to progression %
% (TTP). It also determines protocol that minimizes the normalized AUC  %
% of toxicity.                                                          %
%                                                                       %
% Treatment-related thresholds considered:                              %
% 1) Vary toxicity-related upper bound (turn drug off) and lower bound  %
%    (turn drug on), which is only impactful for protocol #3            %
%    (daily+tox) and protocol #4 (adpative+tox).                        %
% 2) Vary tumor size-related upper bound (turn drug on) and lower bound %
%    (turn drug off), which is only impactful for protocol #2           %
%    (adaptive) and protocol #4 (adpative+tox).                         %
%                                                                       %
% We do not simulate protocols that are infeasible because upper        %
% threshold is <= lower threshold. Those protocols get assigned nan.    %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc; 
[q, ICs] = set_parameters_ICs();

%% Protocol features
tspan = 0:0.1:100;
nrun = max(tspan);           % number of therapy cycles for continuous
tp = 1;                      % therapy period
dose = 1;                    % drug bolus dose
baseline = (ICs(1)+ICs(2));
Nhi = 1.0*baseline;   % adaptive threshold: Rx turns on
Nlo = 0.4*baseline;   % adaptive threshold: Rx turns off 
Nfail = 1.5*baseline; % Rx failure threshold
Thi = 2;                     % Toxicity tolerance: upper bound
Tlo = 1;                     % Toxicity tolerance: lower bound

%% To store output
num_pts = 13; %num_pts = 11; 
tfail_sweep1 = cell(5,1); 
tfail_sweep1_diff = cell(5,1); 
AUC_sweep1 = cell(5,1); 
AUC_sweep1_diff = cell(5,1); 
protocol_string = {};
protocol_string{2} = 'Daily';
protocol_string{3} = 'Adaptive';
protocol_string{4} = 'Daily + Toxicity';
protocol_string{5} = 'Adaptive + Toxicity';
%% Sweep over values of beta and alphaR
beta_default = q.beta;
alphaR_default = q.alphaR;
beta  = linspace(0.5,3.5,num_pts);  % linspace(0.5,4,num_pts); 
alphaR = linspace(0.2*q.alphaS,1*q.alphaS,num_pts); % linspace(0.3*q.alphaS,0.7*q.alphaS,num_pts);

for protocol = 2:5
    fprintf('%s: Sweep in beta-alphaR space\n',protocol_string{protocol}); 
    tfail = zeros(length(beta),length(alphaR)); 
    AUC = zeros(length(beta),length(alphaR)); 
    for i = 1:length(beta)
        q.beta = beta(i);
        for j = 1:length(alphaR)
            q.alphaR = alphaR(j); 
            [Tx,Rx,time,sens,resist,drug,tox] = ...
                solve_model(protocol, ICs, tp, nrun, q, ...
                dose, Nfail, Nlo, Nhi, Thi, Tlo);
            time_to_check = time(find((sens+resist)>Nfail,1));
            if isempty(time_to_check)==1 % it does not progress, set = 150
                time_to_check = max(tspan)*1.5;
            end
            tfail(i,j) = time_to_check;
            AUC(i,j) = trapz(time,tox)/tfail(i,j);
        end
    end
    tfail_sweep1{protocol} = tfail;  
    AUC_sweep1{protocol}   = AUC;
end

tfail_sweep1_toxicity_vs_none = cell(5,1); 
AUC_sweep1_toxicity_vs_none = cell(5,1); 
for protocol = 4:5
    tfail_sweep1_toxicity_vs_none{protocol} = ...
        tfail_sweep1{protocol-2} - tfail_sweep1{protocol}; 
    AUC_sweep1_toxicity_vs_none{protocol} = ...
        AUC_sweep1{protocol-2} - AUC_sweep1{protocol}; 


    % TTP
    max_TTP = max([tfail_sweep1{protocol-2} tfail_sweep1{protocol}],[],'all');
    min_TTP = min([tfail_sweep1{protocol-2} tfail_sweep1{protocol}],[],'all');
    
    % AUC of toxicity plots
    max_AUC = max([AUC_sweep1{protocol-2} AUC_sweep1{protocol}],[],'all');
    min_AUC = min([AUC_sweep1{protocol-2} AUC_sweep1{protocol}],[],'all');

    figure;
    set(gcf, 'Units', 'Normalized','OuterPosition', [0., 0.05, 0.95, 0.95]);
    if protocol == 4
        sgtitle('Daily Treatment','FontSize',20,'FontWeight','bold');
        fname_fig = 'sweep_beta_alphaR_daily';
    elseif protocol == 5
        sgtitle('Adaptive Treatment','FontSize',20,'FontWeight','bold');
        fname_fig = 'sweep_beta_alphaR_adaptive';
    end

    % Time to failure plots
    subplot(2,3,1)
    hHM=heatmap(beta,alphaR,tfail_sweep1{protocol-2}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = '\beta';
    hHM.YLabel = '\alpha_R = \epsilon\alpha_S';
    hHM.Title = 'TTP: No Tox Feedback';
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.ColorLimits = [min_TTP max_TTP];

    subplot(2,3,2)
    hHM=heatmap(beta,alphaR,tfail_sweep1{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = '\beta';
    hHM.YLabel = '\alpha_R = \epsilon\alpha_S';
    hHM.Title = 'TTP: With Tox Feedback';
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.ColorLimits = [min_TTP max_TTP];
        
    subplot(2,3,3)
    hHM=heatmap(beta,alphaR,-1*tfail_sweep1_toxicity_vs_none{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = '\beta';
    hHM.YLabel = '\alpha_R = \epsilon\alpha_S';
    hHM.Title = 'TTP (With Tox Feedback - Without)';    
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;

    subplot(2,3,4)
    hHM=heatmap(beta,alphaR,AUC_sweep1{protocol-2}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = '\beta';
    hHM.YLabel = '\alpha_R = \epsilon\alpha_S';
    hHM.Title = 'Norm Tox AUC: W/o Tox Feedback';
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.CellLabelFormat = '%.1f';
    hHM.ColorLimits = [min_AUC max_AUC];

    subplot(2,3,5)
    hHM=heatmap(beta,alphaR,AUC_sweep1{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = '\beta';
    hHM.YLabel = '\alpha_R = \epsilon\alpha_S';
    hHM.Title = 'Norm Tox AUC: With Tox Feedback';
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.CellLabelFormat = '%.1f';
    hHM.ColorLimits = [min_AUC max_AUC];
        
    subplot(2,3,6)
    hHM=heatmap(beta,alphaR,AUC_sweep1_toxicity_vs_none{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = '\beta';
    hHM.YLabel = '\alpha_R = \epsilon\alpha_S';
    hHM.Title = 'Norm Tox AUC: (Without - With Tox Feedback)';
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.CellLabelFormat = '%.1f';

    saveas(gcf,[fname_fig,'.fig'])
    saveas(gcf,[fname_fig,'.png'])
end

%% Rank protocols by maximum time to progression
protocol_rank1 = zeros(length(beta),length(alphaR)); 
max_TTP1  = zeros(length(beta),length(alphaR)); 
for i = 1:length(beta)
    for j = 1:length(alphaR)
        TTP_ij = zeros(1,4);
        for protocol = 2:5
            TTP_ij(protocol-1) = tfail_sweep1{protocol}(i,j);
        end
        M = max(TTP_ij,[],"all");
        indices = find(TTP_ij == M);
        protocol_rank1(i,j) = classify(indices); 
        max_TTP1(i,j) = M; 
    end
end

%% Map protocol ranks onto minimal colorbar
B = unique(protocol_rank1);
protocol_rank_renumbered1 = zeros(length(beta),length(alphaR)); 
for k = 1:length(B)
    % Replace each value 
    for i = 1:length(beta)
        for j = 1:length(alphaR)
            if protocol_rank1(i,j) == B(k)
                protocol_rank_renumbered1(i,j) = k;
            end
        end
    end    
end
txt = {'1 = Daily','2 = Adaptive','3 = Daily+Toxicity','4 = Adaptive+Toxicity'};

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.6]);
subplot(1,2,1)
imagesc(beta,alphaR,protocol_rank_renumbered1');
set(gca,'YDir','normal') 
colormap(parula)
clim([1, length(B)]);
for row = 1:size(protocol_rank1, 1)
    for col = 1:size(protocol_rank1, 2)
        [font_color,rank_to_print] = assign_label(protocol_rank1(col, row),1);
        text(beta(col), alphaR(row), rank_to_print, 'HorizontalAlignment', 'center', ...
            'Color', font_color, 'FontWeight', 'bold');
    end
end
xlabel('\beta','FontSize',14);
ylabel('\alpha_R = \epsilon\alpha_S','FontSize',14); 
xticks(beta)
yticks(alphaR)
title('Protocol that Maximizes TTP','FontSize',16);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

subplot(1,2,2)
title('Protocol Key')
t = text(0,0.5,txt);
t.FontSize = 14;
set(gca, 'visible', 'off');
fname_fig = 'sweep_beta_alphaR_protocol_ranking';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])


figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.6]);
subplot(1,2,1)
hHM=heatmap(beta,alphaR,max_TTP1');
hHM.YDisplayData=flip(hHM.YDisplayData);
hHM.XLabel = '\beta';
hHM.YLabel = '\alpha_R = \epsilon\alpha_S';
hHM.Title = 'Maximium TTP';
hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
hHM.FontSize = 14;
%hHM.ColorbarVisible = 'off';

subplot(1,2,2)
title('Protocol Key')
t = text(0,0.5,txt);
t.FontSize = 14;
set(gca, 'visible', 'off');
fname_fig = 'sweep_beta_alphaR_protocol_ranking_TTP';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

%% Rank protocols by minimum normalized AUC(toxicity)
protocol_rank2 = zeros(length(beta),length(alphaR)); 
min_AUC  = zeros(length(beta),length(alphaR)); 
for i = 1:length(beta)
    for j = 1:length(alphaR)
        AUC_ij = zeros(1,4);
        for protocol = 2:5
            AUC_ij(protocol-1) = AUC_sweep1{protocol}(i,j);
        end
        M = min(AUC_ij,[],"all");
        indices = find(AUC_ij == M);
        protocol_rank2(i,j) = classify(indices); 
        min_AUC(i,j) = M; 
    end
end

%% Map protocol ranks onto minimal colorbar
B2 = unique(protocol_rank2);
protocol_rank_renumbered2 = zeros(length(beta),length(alphaR)); 
for k = 1:length(B2)
    % Replace each value 
    for i = 1:length(beta)
        for j = 1:length(alphaR)
            if protocol_rank2(i,j) == B2(k)
                protocol_rank_renumbered2(i,j) = k;
            end
        end
    end    
end

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.95, 0.6]);
subplot(1,3,1)
imagesc(beta,alphaR,protocol_rank_renumbered2');
set(gca,'YDir','normal') 
colormap(parula)
clim([1, length(B2)+0.01]);
for row = 1:size(protocol_rank2, 1)
    for col = 1:size(protocol_rank2, 2)
        [font_color,rank_to_print] = assign_label(protocol_rank2(col, row),2);
        text(beta(col), alphaR(row), rank_to_print, 'HorizontalAlignment', 'center', ...
            'Color', font_color, 'FontWeight', 'bold');
    end
end
xlabel('\beta','FontSize',14);
ylabel('\alpha_R = \epsilon\alpha_S','FontSize',14); 
title('Protocol that Minimizes Norm AUC of Tox','FontSize',16);

subplot(1,3,2)
hHM=heatmap(beta,alphaR,min_AUC');
hHM.YDisplayData=flip(hHM.YDisplayData);
hHM.XLabel = '\beta';
hHM.YLabel = '\alpha_R = \epsilon\alpha_S';
hHM.Title = 'Mnimum Norm AUC of Tox';
hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
hHM.FontSize = 14;
%hHM.ColorbarVisible = 'off';

subplot(1,3,3) 
title('Protocol Key')
t = text(0,0.5,txt);
t.FontSize = 14;
set(gca, 'visible', 'off');
fname_fig = 'sweep_beta_alphaR_protocol_ranking_AUC';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

%% Now do protocol sensitivity to thresholds at nominal parameters
% Reset beta and alphaR
q.beta = beta_default;
q.alphaR = alphaR_default;
%% Thi
Thi_vec = linspace(1.5,3.5,num_pts); % must keep above Tlow = 1
tfail_sweep2 = cell(5,1); 
AUC_sweep2 = cell(5,1); 
for protocol = 2:5
    fprintf('%s: Sweep over Thi\n',protocol_string{protocol}); 
    tfail = zeros(size(Thi_vec)); 
    AUC = zeros(size(Thi_vec)); 
    for i = 1:length(Thi_vec)
        [Tx,Rx,time,sens,resist,drug,tox] = solve_model(protocol, ICs,...
            tp, nrun, q, dose, Nfail, Nlo, Nhi, Thi_vec(i), Tlo);
        time_to_check = time(find((sens+resist)>Nfail,1));
        if isempty(time_to_check)==1 % it does not progress, set = 150
            time_to_check = max(tspan)*1.5;
        end
        tfail(i) = time_to_check;
        AUC(i) = trapz(time,tox)/tfail(i);
    end
    tfail_sweep2{protocol} = tfail;  
    AUC_sweep2{protocol}   = AUC;
end

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.9]);
for protocol = 2:1:5
    subplot(2,2,protocol-1)
    % Shading above threshold
    rectangle('Position',[min(Thi_vec),nrun,max(Thi_vec)-min(Thi_vec),...
        nrun*1.55-nrun],'FaceColor',[0.9, 0.9, 0.9],'EdgeColor','none');
    hold on;
    plot(Thi_vec,tfail_sweep2{protocol},'-o','LineWidth',2);
    hold off;
    xlabel('Upper toxicity bound: Tox_o_n','FontSize',14)
    ylabel('Time to progression','FontSize',14);
    title(protocol_string{protocol},'FontSize',18);
    xlim([min(Thi_vec),max(Thi_vec)])
    ylim([0,nrun*1.55]);
end
fname_fig = 'sweep_Tox_off';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

%% Tlo
Tlo_vec = linspace(0.25,1.75,num_pts); % must keep below Thi = 2
tfail_sweep3 = cell(5,1); 
AUC_sweep3 = cell(5,1); 
for protocol = 2:5
    fprintf('%s: Sweep over Tlo\n',protocol_string{protocol}); 
    tfail = zeros(size(Tlo_vec)); 
    AUC = zeros(size(Tlo_vec)); 
    for i = 1:length(Tlo_vec)
        [Tx,Rx,time,sens,resist,drug,tox] = solve_model(protocol, ICs,...
            tp, nrun, q, dose, Nfail, Nlo, Nhi, Thi, Tlo_vec(i));
        time_to_check = time(find((sens+resist)>Nfail,1));
        if isempty(time_to_check)==1 % it does not progress, set = 150
            time_to_check = max(tspan)*1.5;
        end
        tfail(i) = time_to_check;
        AUC(i) = trapz(time,tox)/tfail(i);
    end
    tfail_sweep3{protocol} = tfail;  
    AUC_sweep3{protocol}   = AUC;
end

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.9]);
for protocol = 2:1:5
    subplot(2,2,protocol-1)
    % Shading above threshold
    rectangle('Position',[min(Tlo_vec),nrun,max(Tlo_vec)-min(Tlo_vec),...
        nrun*1.55-nrun],'FaceColor',[0.9, 0.9, 0.9],'EdgeColor','none');
    hold on;
    plot(Tlo_vec,tfail_sweep3{protocol},'-o','LineWidth',2);
    hold off;
    xlabel('Lower toxicity bound: Tox_o_n','FontSize',14)
    ylabel('Time to progression','FontSize',14);
    title(protocol_string{protocol},'FontSize',18);
    xlim([min(Tlo_vec),max(Tlo_vec)])
    ylim([0,nrun*1.55]);
end
fname_fig = 'sweep_Tox_on';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

%% Nhi
Nhi_vec = linspace(25,Nfail,num_pts); % must keep above Nlo = 20 and below Nfail
tfail_sweep4 = cell(5,1); 
AUC_sweep4 = cell(5,1); 
for protocol = 2:5
    fprintf('%s: Sweep over Nhi\n',protocol_string{protocol}); 
    tfail = zeros(size(Nhi_vec)); 
    AUC = zeros(size(Nhi_vec)); 
    for i = 1:length(Nhi_vec)
        [Tx,Rx,time,sens,resist,drug,tox] = solve_model(protocol, ICs,...
            tp, nrun, q, dose, Nfail, Nlo, Nhi_vec(i), Thi, Tlo);
        time_to_check = time(find((sens+resist)>Nfail,1));
        if isempty(time_to_check)==1 % it does not progress, set = 150
            time_to_check = max(tspan)*1.5;
        end
        tfail(i) = time_to_check;
        AUC(i) = trapz(time,tox)/tfail(i);
    end
    tfail_sweep4{protocol} = tfail;  
    AUC_sweep4{protocol}   = AUC;
end

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.9]);
for protocol = 2:1:5
    subplot(2,2,protocol-1)
    % Shading above threshold
    rectangle('Position',[min(Nhi_vec),nrun,...
        max(Nhi_vec)-min(Nhi_vec),...
        nrun*1.55-nrun],'FaceColor',[0.9, 0.9, 0.9],'EdgeColor','none');
    hold on;
    plot(Nhi_vec,tfail_sweep4{protocol},'-o','LineWidth',2);
    hold off;
    xlabel('Upper size bound: Rx_o_n','FontSize',14)
    ylabel('Time to progression','FontSize',14);
    title(protocol_string{protocol},'FontSize',18);
    xlim([min(Nhi_vec),max(Nhi_vec)])
    ylim([0,nrun*1.55]);
end
fname_fig = 'sweep_Rx_on';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

%% Nlo
Nlo_vec = linspace(5,45,num_pts); % must keep below Nhi = 50
tfail_sweep5 = cell(5,1); 
AUC_sweep5 = cell(5,1); 
for protocol = 2:5
    fprintf('%s: Sweep over Nlo\n',protocol_string{protocol}); 
    tfail = zeros(size(Nlo_vec)); 
    AUC = zeros(size(Nlo_vec)); 
    for i = 1:length(Nlo_vec)
        [Tx,Rx,time,sens,resist,drug,tox] = solve_model(protocol, ICs,...
            tp, nrun, q, dose, Nfail, Nlo_vec(i), Nhi, Thi, Tlo);
        time_to_check = time(find((sens+resist)>Nfail,1));
        if isempty(time_to_check)==1 % it does not progress, set = 150
            time_to_check = max(tspan)*1.5;
        end
        tfail(i) = time_to_check;
        AUC(i) = trapz(time,tox)/tfail(i);
    end
    tfail_sweep5{protocol} = tfail;  
    AUC_sweep5{protocol}   = AUC;
end

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.9]);
for protocol = 2:1:5
    subplot(2,2,protocol-1)
    % Shading above threshold
    rectangle('Position',[min(Nlo_vec),nrun,...
        max(Nlo_vec)-min(Nlo_vec),...
        nrun*1.55-nrun],'FaceColor',[0.9, 0.9, 0.9],'EdgeColor','none');
    hold on;
    plot(Nlo_vec,tfail_sweep5{protocol},'-o','LineWidth',2);
    hold off;
    xlabel('Lower size bound: Rx_o_f_f','FontSize',14)
    ylabel('Time to progression','FontSize',14);
    title(protocol_string{protocol},'FontSize',18);
    xlim([min(Nlo_vec),max(Nlo_vec)])
    ylim([0,nrun*1.55]);
end
fname_fig = 'sweep_Rx_off';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])
 

%% Sweep over values of Tlo and Thi
Tlo_vec = linspace(0.25,2.25,num_pts); 
Thi_vec = linspace(0.5,4,num_pts); 
tfail_sweep6 = cell(5,1); 
AUC_sweep6 = cell(5,1); 
for protocol = 2:5
    fprintf('%s: Sweep in Tlo-Thi space\n',protocol_string{protocol}); 
    tfail = zeros(length(Tlo_vec),length(Thi_vec)); 
    AUC = zeros(length(Tlo_vec),length(Thi_vec)); 
    for i = 1:length(Tlo_vec)
        for j = 1:length(Thi_vec)
            if Thi_vec(j)<=Tlo_vec(i) % not allowable protocol
                tfail(i,j) = nan;
                AUC(i,j) = nan;
            else
                [Tx,Rx,time,sens,resist,drug,tox] = ...
                    solve_model(protocol, ICs, tp, nrun, q, ...
                    dose, Nfail, Nlo, Nhi, Thi_vec(j), Tlo_vec(i));
                time_to_check = time(find((sens+resist)>Nfail,1));
                if isempty(time_to_check)==1 % it does not progress, set = 150
                    time_to_check = max(tspan)*1.5;
                end
                tfail(i,j) = time_to_check;
                AUC(i,j) = trapz(time,tox)/tfail(i,j);
            end
        end
    end
    tfail_sweep6{protocol} = tfail;  
    AUC_sweep6{protocol}   = AUC;
end

figure;
sgtitle('Toxicity Thresholds','FontSize',18,'FontWeight','bold');
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.55, 0.95]);
for protocol = 2:5
    % Time to failure plots
    subplot(2,2,protocol-1)
    hHM=heatmap(Tlo_vec,Thi_vec,tfail_sweep6{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = 'Lower toxicity bound: Tox_o_n';
    hHM.YLabel = 'Upper toxicity bound: Tox_o_f_f';
    hHM.Title = protocol_string{protocol};
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.ColorLimits = [0 150];
end
fname_fig = 'sweep_ToxOn_ToxOff_TTP_all';    
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

figure;
sgtitle('Toxicity Thresholds','FontSize',18,'FontWeight','bold');
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.7]);
for protocol = 4:5
    % Time to failure plots
    subplot(1,2,protocol-3)
    hHM=heatmap(Tlo_vec,Thi_vec,tfail_sweep6{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = 'Lower toxicity bound: Tox_o_n';
    hHM.YLabel = 'Upper toxicity bound: Tox_o_f_f';
    title_append = num2str(tfail_sweep6{protocol-2}(1,1));
    hHM.Title = [protocol_string{protocol} ': TTP(no tox) = ' title_append];
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.ColorLimits = [0 150];
end
fname_fig = 'sweep_ToxOn_ToxOff_TTP_toxOnly';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

%% Sweep over values of Nlo and Nhi
Nhi_vec = linspace(25,Nfail,num_pts); % must keep above Nlo = 20 and below Nfail
Nlo_vec = linspace(5,45,num_pts); 
tfail_sweep7 = cell(5,1); 
AUC_sweep7 = cell(5,1); 
for protocol = 2:5
    fprintf('%s: Sweep in Nlo-Nhi space\n',protocol_string{protocol}); 
    tfail = zeros(length(Nlo_vec),length(Nhi_vec)); 
    AUC = zeros(length(Nlo_vec),length(Nhi_vec)); 
    for i = 1:length(Nhi_vec)
        for j = 1:length(Nlo_vec)
            if Nhi_vec(i)<=Nlo_vec(j) % not allowable protocol
                tfail(i,j) = nan;
                AUC(i,j) = nan;
            else
                [Tx,Rx,time,sens,resist,drug,tox] = ...
                    solve_model(protocol, ICs, tp, nrun, q, ...
                    dose, Nfail, Nlo_vec(j), Nhi_vec(i), Thi, Tlo);
                time_to_check = time(find((sens+resist)>Nfail,1));
                if isempty(time_to_check)==1 % it does not progress, set = 150
                    time_to_check = max(tspan)*1.5;
                end
                tfail(i,j) = time_to_check;
                AUC(i,j) = trapz(time,tox)/tfail(i,j);
            end
        end
    end
    tfail_sweep7{protocol} = tfail;  
    AUC_sweep7{protocol}   = AUC;
end

figure;
sgtitle('Adaptive Thresholds','FontSize',18,'FontWeight','bold');
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.55, 0.95]);
for protocol = 2:5
    % Time to failure plots
    subplot(2,2,protocol-1)
    hHM=heatmap(Nhi_vec,Nlo_vec,tfail_sweep7{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = 'Upper size bound: Rx_o_n';
    hHM.YLabel = 'Lower size bound: Rx_o_f_f';
    hHM.Title = protocol_string{protocol};
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.ColorLimits = [0 150];
end
fname_fig = 'sweep_RxOn_RxOff_TTP_all';    
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

figure;
sgtitle('Adaptive Thresholds','FontSize',18,'FontWeight','bold');
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.7]);
for protocol = 3:2:5
    % Time to failure plots
    if protocol == 3
        subplot(1,2,1)
        title_append = num2str(tfail_sweep7{2}(1,1)); % compare to daily
    elseif protocol == 5
        subplot(1,2,2)
        title_append = num2str(tfail_sweep7{4}(1,1)); % compare to daily+tox
    end
    hHM=heatmap(Nhi_vec,Nlo_vec,tfail_sweep7{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = 'Upper size bound: Rx_o_n';
    hHM.YLabel = 'Lower size bound: Rx_o_f_f';
    if protocol == 3
        hHM.Title = [protocol_string{protocol} ': TTP(daily) = ' title_append];
    elseif protocol == 5
        hHM.Title = [protocol_string{protocol} ': TTP(daily+tox) = ' title_append];
    end
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.ColorLimits = [0 150];
end
fname_fig = 'sweep_RxOn_RxOff_TTP_adaptOnly';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])
 
save output_sweeps.mat beta alphaR protocol_string tfail_sweep1...
    tfail_sweep1_toxicity_vs_none protocol_rank1 AUC_sweep1 ...
    AUC_sweep1_toxicity_vs_none protocol_rank2 tfail_sweep2 ...
    tfail_sweep3 tfail_sweep4 tfail_sweep5 tfail_sweep6 tfail_sweep7

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setup function
function [q, ICs] = set_parameters_ICs()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%     Tumor Growth/Death Params     %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    q.alphaS = 1;           % per unit time, sensitive cell growth rate (fixed)
    q.epsilon = 0.4; 
    q.alphaR =  q.epsilon*q.alphaS;   % per unit time, resistant cell growth rate
                            % alphaR > 1 => R grows faster than S
                            % alphaR < 1 => R grows slower than S
                            % alphaR = 1 => R and S have same proliferation
                            % rate
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


%% Tumor model equations
function dxdt = tumor_model(t, x, q)
    S = x(1);  % Sensitive cells
    R = x(2);  % Resistant cells
    D = x(3);  % Drug concentration
    T = x(4);  % Toxicity level

    % Tumor and drug dynamics
    dS = q.alphaS*S*(1 - (S + R)/q.K) - q.delta*S*D;
    dR = q.alphaR*R*(1 - (q.beta*S + R)/q.K); 
    dD = - q.lambda*D;
    dT = q.mu*D - q.gamma*T;
    
    dxdt = [dS; dR; dD; dT];
end

%% Main dosing function with feedbacks
function [Tx,Rx,time,sens,resist,drug,tox] = solve_model(protocol, ...
    ICs, tp, nrun, q, dose, Nfail, Nlo, Nhi, Thi, Tlo)
    time = []; drug = []; tox = []; 
    sens = []; resist = []; 
    Rx = zeros(1,nrun);
    Tx = zeros(1,nrun);
    ICs(3) = dose; 
  
    on = 1;   
    for i = 1:1:nrun
        trun = (i-1)*tp:0.1:i*tp;
        Rx(i) = on*150;         % to draw gray background for Rx
        Tx(i) = trun(end)-1;    % to draw gray background for Rx

        % Solve ODEs
        [t,F] = ode23s(@(t, x) tumor_model(t, x, q), trun, ICs);

        % update time and cell arrays
        time = vertcat(time,t(2:end));
        drug = vertcat(drug,F(2:end,3));
        sens = vertcat(sens,F(2:end,1));
        resist = vertcat(resist,F(2:end,2));
        tox = vertcat(tox,F(2:end,4));
    
        on_new = on; % by default, we don't change the strategy
        if protocol == 2 % daily, no toxicity feedback
            % Nothing to check, as no tumor size or toxicity feedback
        elseif protocol == 3 % adaptive, no toxicity feedback
            % check to see if tumor size condition for off cycling met
            %N = F(:, 1) + F(:, 2);
            Nend = F(end, 1) + F(end, 2);
            if on == 1 && (Nend < Nlo) %any(N < Nlo)
                on_new = 0;
                % fprintf('\tTurning drug OFF at t = %f since tumor size = %f < %f\n',...
                %     max(trun),Nend,Nlo);
            end
            if on == 0 && (Nend > Nhi) %any(N > Nhi)
                on_new = 1;
                % fprintf('\tTurning drug ON at t = %f since tumor size = %f > %f\n',...
                %     max(trun),Nend,Nhi);
            end
        elseif protocol == 4 % daily + toxicity feedback
            % check to see if toxicity condition for off cycling 
            tox_end = F(end,4); 
            if (on == 1) && (tox_end > Thi) 
                on_new = 0;
                % fprintf('\tTurning drug OFF at t = %f since end(tox) = %f > %f\n',...
                %     max(trun),tox_end,Thi);
            end
            if (on == 0) && (tox_end < Tlo)
                on_new = 1;
                % fprintf('\tTurning drug ON at t = %f since end(tox) = %f < %f\n',...
                %     max(trun),tox_end,Thi);
            end
        elseif protocol == 5 % adaptive + toxicity feedback
            Nend = F(end, 1) + F(end, 2);
            tox_end = F(end,4); 
            if (on == 1) && ((Nend < Nlo) || (tox_end > Thi))
                on_new = 0;
                % fprintf('\tTurning drug OFF at t = %f for one of two reasons:\n',max(trun));
                % fprintf('\t\tEither tumor size = %f < %f,\n',Nend,Nlo);
                % fprintf('\t\tOr, end(tox) = %f > %f\n',tox_end,Thi);
            end
            if (on == 0) && (Nend > Nhi) && (tox_end < Tlo)
                on_new = 1;
                % fprintf('\tTurning drug ON at t = %f since both:\n',max(trun));
                % fprintf('\t\tTumor size = %f > %f,\n',Nend,Nhi);
                % fprintf('\t\tAnd, end(tox) = %f < %f\n',tox_end,Tlo);
            end
        end
        on = on_new;

        
        ICs = [F(end,1), F(end,2), F(end,3)+on*dose, F(end,4)]; 
        clear t F 
        if max(sens+resist)>Nfail
            break
        end      
    end  
end

function protocol_rank = classify(indices)
    %display(indices)
    if length(indices) == 1
        if indices(1) == 1
            protocol_rank = 1;
        elseif indices(1) == 2
            protocol_rank = 2;
        elseif indices(1) == 3
            protocol_rank = 3;
        elseif indices(1) == 4
            protocol_rank = 4;
        end
    elseif length(indices) == 2
        if indices(1) == 1 && indices(2) == 2
            protocol_rank = 5;
        elseif indices(1) == 1 && indices(2) == 3
            protocol_rank = 6;
        elseif indices(1) == 1 && indices(2) == 4
            protocol_rank = 7;
        elseif indices(1) == 2 && indices(2) == 3
            protocol_rank = 8;
        elseif indices(1) == 2 && indices(2) == 4
            protocol_rank = 9;
        elseif indices(1) == 3 && indices(2) == 4
            protocol_rank = 10; 
        end
    elseif length(indices) == 3
        if indices(1) == 1 && indices(2) == 2 && indices(3) == 3
            protocol_rank = 11;
        elseif indices(1) == 1 && indices(2) == 2 && indices(3) == 4
            protocol_rank = 12;
        elseif indices(1) == 1 && indices(2) == 3 && indices(3) == 4
            protocol_rank = 13;
        elseif indices(1) == 2 && indices(2) == 3 && indices(3) == 4
            protocol_rank = 14;
        end
    elseif length(indices) == 4
        protocol_rank = 15; % all protocols give same TTP
    else
        fprintf('length(indices) = %d, but should be between 1 and 4\n',...
            length(indices)); 
        STOP
    end
    %fprintf('\tProtocol is classified as rank %d\n',protocol_rank);        
end

function [font_color,rank_to_print] = assign_label(protocol_rank_ij,flag)
        font_color='w';
        if protocol_rank_ij <= 4
            rank_to_print = num2str(protocol_rank_ij);
        elseif protocol_rank_ij == 5
            rank_to_print = '1/2';
        elseif protocol_rank_ij == 6
            rank_to_print = '1/3';
        elseif protocol_rank_ij == 7
            rank_to_print = '1/4';
        elseif protocol_rank_ij == 8
            rank_to_print = '2/3';
        elseif protocol_rank_ij == 9
            rank_to_print = '2/4';
        elseif protocol_rank_ij== 10
            rank_to_print = '3/4';
            if flag == 2
                font_color = 'k';
            end
        elseif protocol_rank_ij == 11
            rank_to_print = '1/2/3';
        elseif protocol_rank_ij == 12
            rank_to_print = '1/2/4';
        elseif protocol_rank_ij == 13
            rank_to_print = '1/3/4';
        elseif protocol_rank_ij == 14
            rank_to_print = '2/3/4';
            if flag == 1
                font_color = 'k';
            end
        elseif protocol_rank_ij == 15
            rank_to_print = '1/2/3/4';
        end
end