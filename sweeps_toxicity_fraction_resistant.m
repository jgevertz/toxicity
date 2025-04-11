%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Quantifying the impact of model parameters on the REASON treatment    %
% fails. That is, is failure a result of sensitive cells not being 		%
% adequately targeted by the drug, or a result of resistance? 			%
% Updated: 2/24/2025.                                                   %
%                                                                       %
% For several values of the sesnsitive cell growth rate (alphaS), sweep %
% over the following parameters:  										%
% - Competition parameter, beta                                         %
% - Drug-induced kill rate of sensitive cells, delta                    %
%                                                                       %
% The role of the toxicity decay rate (gamma) is then explored at one 	%
% parameterization for which treatment fails because of sensitive cells %
% and at another for which treatment fails because of resistant cells. 	%
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc; 
[q, ICs] = set_parameters_ICs();
%q.alphaS = 1.45; % study impact of sensitive cell growth rate
%q.alphaR = q.epsilon*q.alphaS; 

%% Protocol features
tspan = 0:0.1:100;
nrun = max(tspan);           % number of therapy cycles for continuous
tp = 1;                      % therapy period
dose = 1;                    % drug bolus dose
baseline = (ICs(1)+ICs(2));
Nhi = 1.0*baseline;   % adaptive threshold: Rx turns on
Nlo = 0.4*baseline;   % adaptive threshold: Rx turns off 
Nfail = 1.5 * baseline; % Rx failure threshold
Thi = 2;                     % Toxicity tolerance: upper bound
Tlo = 1;                     % Toxicity tolerance: lower bound

%% To store output
num_pts = 13; 
tfail_sweep1 = cell(5,1); 
tfail_sweep1_diff = cell(5,1); 
AUC_sweep1 = cell(5,1); 
AUC_sweep1_diff = cell(5,1); 
Rfrac_fail_sweep1 = cell(5,1);
protocol_string = {};
protocol_string{2} = 'Daily';
protocol_string{3} = 'Adaptive';
protocol_string{4} = 'Daily + Toxicity';
protocol_string{5} = 'Adaptive + Toxicity';
%% Sweep over values of beta and delta
beta_default = q.beta;
delta_default = q.delta;
beta  = linspace(0.5,3.5,num_pts);
delta = linspace(0.5,1.5,num_pts);
for protocol = 2:5
    fprintf('%s: Sweep in beta-delta space\n',protocol_string{protocol}); 
    tfail = zeros(length(beta),length(delta)); 
    AUC = zeros(length(beta),length(delta)); 
    Rfrac_fail = zeros(length(beta),length(delta)); 
    for i = 1:length(beta)
        q.beta = beta(i);
        for j = 1:length(delta)
            q.delta = delta(j); 
            [Tx,Rx,time,sens,resist,drug,tox] = ...
                solve_model(protocol, ICs, tp, nrun, q, ...
                dose, Nfail, Nlo, Nhi, Thi, Tlo);
         
            time_to_check = time(find((sens+resist)>Nfail,1));
            if isempty(time_to_check)==1 % it does not progress, set = 150
                time_to_check = max(tspan)*1.5;
                Rfrac_fail(i,j) = nan; % does not fail 
                %fprintf('\tbeta = %f, delta = %f does NOT progress so Rfrac = %f\n',...
                %    beta(i),delta(j),Rfrac_fail(i,j));
            else 
                time_idx = find((sens+resist)>Nfail,1); % time index at failure
                denom = sens(time_idx)+resist(time_idx);
                Rfrac_fail(i,j) = resist(time_idx)/denom; 
                %fprintf('\tbeta = %f, delta = %f progresses with TTP = %f so Rfrac = %f\n',...
                %    beta(i),delta(j),time_to_check,Rfrac_fail(i,j));
            end
            tfail(i,j) = time_to_check;
            AUC(i,j) = trapz(time,tox)/tfail(i,j);
        end
    end
    tfail_sweep1{protocol} = tfail;  
    AUC_sweep1{protocol}   = AUC;
    Rfrac_fail_sweep1{protocol}   = Rfrac_fail;
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
    set(gcf, 'Units', 'Normalized','OuterPosition', [0., 0.05, 1, 0.95]);
    if protocol == 4
        sgtitle('Daily Treatment','FontSize',20,'FontWeight','bold');
        fname_fig = ['sweep_beta_delta_daily_alphaS_' num2str(q.alphaS) ...
            '_eps_' num2str(q.epsilon)];
    elseif protocol == 5
        sgtitle('Adaptive Treatment','FontSize',20,'FontWeight','bold');
        fname_fig = ['sweep_beta_delta_adaptive_' num2str(q.alphaS)...
            '_eps_' num2str(q.epsilon)];
    end

    % Time to failure plots
    subplot(2,3,1)
    hHM=heatmap(beta,delta,tfail_sweep1{protocol-2}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = '\beta';
    hHM.YLabel = '\delta';
    hHM.Title = 'TTP: No Tox Feedback';
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.ColorLimits = [min_TTP max_TTP];

    subplot(2,3,2)
    hHM=heatmap(beta,delta,tfail_sweep1{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = '\beta';
    hHM.YLabel = '\delta';
    hHM.Title = 'TTP: With Tox Feedback';
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.ColorLimits = [min_TTP max_TTP];
        
    subplot(2,3,3)
    hHM=heatmap(beta,delta,-1*tfail_sweep1_toxicity_vs_none{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = '\beta';
    hHM.YLabel = '\delta';
    hHM.Title = 'TTP (With Tox Feedback - Without)';    
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;

    subplot(2,3,4)
    hHM=heatmap(beta,delta,AUC_sweep1{protocol-2}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = '\beta';
    hHM.YLabel = '\delta';
    hHM.Title = 'Norm Tox AUC: W/o Tox Feedback';
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.CellLabelFormat = '%.1f';
    hHM.ColorLimits = [min_AUC max_AUC];

    subplot(2,3,5)
    hHM=heatmap(beta,delta,AUC_sweep1{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = '\beta';
    hHM.YLabel = '\delta';
    hHM.Title = 'Norm Tox AUC: With Tox Feedback';
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.CellLabelFormat = '%.1f';
    hHM.ColorLimits = [min_AUC max_AUC];
        
    subplot(2,3,6)
    hHM=heatmap(beta,delta,AUC_sweep1_toxicity_vs_none{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.XLabel = '\beta';
    hHM.YLabel = '\delta';
    hHM.Title = 'Norm Tox AUC: (Without - With Tox Feedback)';
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.CellLabelFormat = '%.1f';

    saveas(gcf,[fname_fig,'.fig'])
    saveas(gcf,[fname_fig,'.png'])
end

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0., 0.05, 0.75, 0.95]);
sgtitle('Fraction Resistant Cells at TTP','FontSize',20,'FontWeight','bold');
for protocol = 2:5    
    % Time to failure plots
    subplot(2,2,protocol-1)
    hHM=heatmap(beta,delta,Rfrac_fail_sweep1{protocol}');
    hHM.YDisplayData=flip(hHM.YDisplayData);
    hHM.CellLabelFormat = '%.3f';
    hHM.XLabel = '\beta';
    hHM.YLabel = '\delta';
    if protocol == 2
        hHM.Title = 'Daily: No Tox Feedback';
    elseif protocol == 3
        hHM.Title = 'Adaptive: No Tox Feedback';
    elseif protocol == 4
        hHM.Title = 'Daily: With Tox Feedback';
    elseif protocol == 5
        hHM.Title = 'Adaptive: With Tox Feedback';
    end
    hHM.XDisplayLabels = compose('%.2f',str2double(hHM.XDisplayLabels));
    hHM.YDisplayLabels = compose('%.2f',str2double(hHM.YDisplayLabels));
    hHM.FontSize = 14;
    hHM.ColorLimits = [0 1];
end
fname_fig = ['sweep_beta_delta_resistant_fraction_' num2str(q.alphaS)...
            '_eps_' num2str(q.epsilon)];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

%% How does gamma impact a patient that fails because of resistant cells?
q.beta = 2.25;
q.delta = 4/3;
gamma = linspace(0.2,0.8,3*num_pts);
tfail_sweep2 = cell(5,1); 
AUC_sweep2 = cell(5,1); 
Rfrac_fail_sweep2 = cell(5,1);
for protocol = 2:5
    tfail = zeros(size(gamma)); 
    AUC = zeros(size(gamma)); 
    Rfrac_fail = zeros(size(gamma)); 
    for i = 1:length(gamma)
        q.gamma = gamma(i);
        [Tx,Rx,time,sens,resist,drug,tox] = ...
            solve_model(protocol, ICs, tp, nrun, q, dose, Nfail, Nlo, Nhi, Thi, Tlo);
     
        time_to_check = time(find((sens+resist)>Nfail,1));
        if isempty(time_to_check)==1 % it does not progress, set = 150
            time_to_check = max(tspan)*1.5;
            Rfrac_fail(i) = nan; % does not fail 
            %fprintf('\tbeta = %f, delta = %f does NOT progress so Rfrac = %f\n',...
            %    beta(i),delta(j),Rfrac_fail(i,j));
        else 
            time_idx = find((sens+resist)>Nfail,1); % time index at failure
            denom = sens(time_idx)+resist(time_idx);
            Rfrac_fail(i) = resist(time_idx)/denom; 
            %fprintf('\tbeta = %f, delta = %f progresses with TTP = %f so Rfrac = %f\n',...
            %    beta(i),delta(j),time_to_check,Rfrac_fail(i,j));
        end
        tfail(i) = time_to_check;
        AUC(i) = trapz(time,tox)/tfail(i);
    end
    tfail_sweep2{protocol} = tfail;  
    AUC_sweep2{protocol}   = AUC;
    Rfrac_fail_sweep2{protocol}   = Rfrac_fail;
end

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0., 0.05, 0.55, 0.85]);
sgtitle(['Failure due to Resistant Cells: \delta = ' num2str(q.delta)],...
    'FontSize',20,'FontWeight','bold');
for protocol = 2:5    
    % Time to failure plots
    subplot(2,2,protocol-1)
    plot(gamma,Rfrac_fail_sweep2{protocol},'o-','LineWidth',2);
    if protocol == 2
        title('Daily: No Tox Feedback','FontSize',16);
    elseif protocol == 3
        title('Adaptive: No Tox Feedback','FontSize',16);
    elseif protocol == 4
        title('Daily: With Tox Feedback','FontSize',16);
    elseif protocol == 5
        title('Adaptive: With Tox Feedback','FontSize',16);
    end
    xlabel('\gamma','FontSize',14);
    ylabel('Fraction Resistant Cells at TTP','FontSize',16)
    ylim([0,1])
end
fname_fig = ['sweep_gamma_delta_eq' num2str(q.delta)];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

%% How does gamma impact a patient that fails because of sensitive cells?
q.delta = 2/3;
tfail_sweep3 = cell(5,1); 
AUC_sweep3 = cell(5,1); 
Rfrac_fail_sweep3 = cell(5,1);
for protocol = 2:5
    tfail = zeros(size(gamma)); 
    AUC = zeros(size(gamma)); 
    Rfrac_fail = zeros(size(gamma)); 
    for i = 1:length(gamma)
        q.gamma = gamma(i);
        [Tx,Rx,time,sens,resist,drug,tox] = ...
            solve_model(protocol, ICs, tp, nrun, q, dose, Nfail, Nlo, Nhi, Thi, Tlo);
     
        time_to_check = time(find((sens+resist)>Nfail,1));
        if isempty(time_to_check)==1 % it does not progress, set = 150
            time_to_check = max(tspan)*1.5;
            Rfrac_fail(i) = nan; % does not fail 
            %fprintf('\tbeta = %f, delta = %f does NOT progress so Rfrac = %f\n',...
            %    beta(i),delta(j),Rfrac_fail(i,j));
        else 
            time_idx = find((sens+resist)>Nfail,1); % time index at failure
            denom = sens(time_idx)+resist(time_idx);
            Rfrac_fail(i) = resist(time_idx)/denom; 
            %fprintf('\tbeta = %f, delta = %f progresses with TTP = %f so Rfrac = %f\n',...
            %    beta(i),delta(j),time_to_check,Rfrac_fail(i,j));
        end
        tfail(i) = time_to_check;
        AUC(i) = trapz(time,tox)/tfail(i);
    end
    tfail_sweep3{protocol} = tfail;  
    AUC_sweep3{protocol}   = AUC;
    Rfrac_fail_sweep3{protocol}   = Rfrac_fail;
end

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0., 0.05, 0.55, 0.85]);
sgtitle(['Failure due to Sensitive Cells: \delta = ' num2str(q.delta)],...
    'FontSize',20,'FontWeight','bold');
for protocol = 2:5    
    % Time to failure plots
    subplot(2,2,protocol-1)
    plot(gamma,Rfrac_fail_sweep3{protocol},'o-','LineWidth',2);
    if protocol == 2
        title('Daily: No Tox Feedback','FontSize',16);
    elseif protocol == 3
        title('Adaptive: No Tox Feedback','FontSize',16);
    elseif protocol == 4
        title('Daily: With Tox Feedback','FontSize',16);
    elseif protocol == 5
        title('Adaptive: With Tox Feedback','FontSize',16);
    end
    xlabel('\gamma','FontSize',14);
    ylabel('Fraction Resistant Cells at TTP','FontSize',16)
    ylim([0,1])
end
fname_fig = ['sweep_gamma_delta_eq' num2str(q.delta)];
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

figure; 
set(gcf, 'Units', 'Normalized','OuterPosition', [0., 0.05, 0.55, 0.85]);
% First plot failure due to resistant cells
plot_num = 1;
for protocol = 4:5
    subplot(2,2,plot_num)
    % Color area where fail due to sensitive cells grey
    rectangle('Position',[min(gamma),0,gamma(19)-min(gamma),1],'FaceColor',...
        [0.9, 0.9, 0.9],'EdgeColor','none');
    hold on;
    % Color area where we don't fail light green
    rectangle('Position',[gamma(19),0,gamma(25)-gamma(19),1],'FaceColor',...
        [212/255, 1, 212/255, 0.5],'EdgeColor','none');
    % Color area where fail due to resistant cells pink
    rectangle('Position',[gamma(25),0,gamma(end)-gamma(25),1],'FaceColor',...
        [1, 0.93, 1,0.5],'EdgeColor','none');
    plot(gamma,Rfrac_fail_sweep3{protocol},'o-','LineWidth',2);
    hold off;
    if protocol == 4
        title('Daily: With Tox Feedback: \delta = 2/3','FontSize',16);
    elseif protocol == 5
        title('Adaptive: With Tox Feedback: \delta = 2/3','FontSize',16);
    end
    plot_num = plot_num + 1; 
    xlabel('\gamma','FontSize',14);
    ylabel('Fraction Resistant Cells at TTP','FontSize',16)
    xlim([gamma(1) gamma(end)])
    ylim([0,1])
end
% Next plot failure due to sensitive cells
gamma_step = gamma(end)-gamma(end-1);
for protocol = 4:5
    subplot(2,2,plot_num)
    % Color area where fail due to sensitive cells grey
    rectangle('Position',[min(gamma),0,gamma(7)-min(gamma)+0.5*gamma_step,1],...
        'FaceColor',[0.9, 0.9, 0.9],'EdgeColor','none');
    % Color area where fail due to resistant cells pink
    rectangle('Position',[gamma(8)-0.5*gamma_step,0,gamma(end)-gamma(8)+0.5*gamma_step,1],...
        'FaceColor',[1, 0.93, 1,0.5],'EdgeColor','none');
    hold on;
    plot(gamma,Rfrac_fail_sweep2{protocol},'o-','LineWidth',2);
    hold off;
    if protocol == 4
        title('Daily: With Tox Feedback: \delta = 4/3','FontSize',16);
    elseif protocol == 5
        title('Adaptive: With Tox Feedback: \delta = 4/3','FontSize',16);
    end
    plot_num = plot_num + 1; 
    xlabel('\gamma','FontSize',14);
    ylabel('Fraction Resistant Cells at TTP','FontSize',16)
    xlim([gamma(1) gamma(end)])
    ylim([0,1])
end
fname_fig = 'sweep_gamma';
saveas(gcf,[fname_fig,'.fig'])
saveas(gcf,[fname_fig,'.png'])

save output_sweeps_resistFrac.mat beta delta tfail_sweep1 ...
    tfail_sweep1_toxicity_vs_none AUC_sweep1 AUC_sweep1_toxicity_vs_none ...
    Rfrac_fail_sweep1 Rfrac_fail_sweep2 Rfrac_fail_sweep3

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