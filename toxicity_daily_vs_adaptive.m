%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Toxicity model simulations at nominal parameterization.               %
% Updated: 2/24/2025                                                    %
% - Toxicity is checked only at the end of window (at dosing decision)  %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc; 
[q, ICs] = set_parameters_ICs();

%% Protocol features
tspan = 0:0.1:100;
nrun = max(tspan);           % number of therapy cycles for continuous
tp = 1;                      % therapy period
dose = 1;                    % drug bolus dose
baseline = (ICs(1) + ICs(2));
Nhi = 1. * baseline;        % adaptive threshold: Rx turns on
Nlo = 0.4 * baseline;        % adaptive threshold: Rx turns off: MUST have Nlo<Nhi
Nfail = 1.5 * baseline;      % Rx failure threshold
Thi = 2;                     % Toxicity tolerance: upper bound: MUST have Thi>Tlo
Tlo = 1;                     % Toxicity tolerance: lower bound 

%% To store output
time_all   = cell(5,1); % Protocol 1 = control, 2 = daily, 3 = adaptive
                        % 4 = daily+toxicity, 5 = adaptive+toxicity
sens_all   = cell(5,1); 
resist_all = cell(5,1); 
drug_all   = cell(5,1); 
tox_all    = cell(5,1); 
tfail_all  = zeros(5,1); 
Tx         = cell(5,1); 
Rx         = cell(5,1); 

%% Control: no drug
[time_all{1},tum_control] = ode23s(@(t, x) tumor_model(t, x, q), tspan, ICs);
sens_all{1}   = tum_control(:,1);
resist_all{1} = tum_control(:,2); 
drug_all{1}   = tum_control(:,3); 
tox_all{1}    = tum_control(:,4); 
tfail_all(1) = time_all{1}(find((sens_all{1}+resist_all{1})>Nfail,1));
% figure; % Plot control
% hold on
% set(gca,'LineWidth',1.25,'FontSize',16,'FontWeight','normal','FontName','Helvetica')
% h1 = plot(time_all{1},sens_all{1},'LineWidth',2,'Color',"#0072BD");
% h2 = plot(time_all{1},resist_all{1},'LineWidth',2,'Color',"#D95319");
% h3 = plot(time_all{1},sens_all{1}+resist_all{1},'--','LineWidth',2,'Color',"#A2142F");
% xlabel('Time','FontSize',16)
% ylabel('Tumor Size','FontSize',16)
% yline(Nfail,'k--','Rx Failure','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','middle')
% yline(Nhi,'k--','Rx On','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','middle')
% yline(Nlo,'k--','Rx Off','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','middle')
% legend([h1,h2,h3],' Sensitive Cells',' Resistant Cells',' Total Cells');
% hold off;
% ylim([0,inf])
% legend boxoff
% title('Control (no treatment)','FontSize',16);


%% Sweep over protocols
% Daily treatment with no toxicity feedback (protocol_switch = 2)
% Adaptive treatment with no toxicity feedback (protocol_switch = 3)
% Daily treatment with toxicity feedback (protocol_switch = 4)
% Adaptive treatment with toxicity feedback (protocol_switch = 5)
for i = 2:5 
    [Tx{i},Rx{i},time_all{i},sens_all{i},resist_all{i},drug_all{i},tox_all{i}] = ...
        solve_model(i, ICs, tp, nrun, q, dose, Nfail, Nlo, Nhi, Thi, Tlo);
    time_to_check = time_all{i}(find((sens_all{i}+resist_all{i})>Nfail,1));
    if isempty(time_to_check)==1 % it does not progress, set = 150
        time_to_check = max(tspan)*1.5;
    end
    tfail_all(i) = time_to_check;
    if i == 2
        protocol_string = 'daily';
    elseif i == 3
        protocol_string = 'adaptive';
    elseif i == 4
        protocol_string = 'daily+toxicity';
    elseif i == 5
        protocol_string = 'adaptive+toxicity';
    end
    fprintf('\t\tFor %s protocol, time until treatment failure = %f\n',...
        protocol_string,tfail_all(i));
end

%% Plot results
for i = 2:5 
    [x,y] = stairs(Tx{i},Rx{i}); 
    figure; 
    set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.8, 0.7]);
    if i == 2
        sgtitle('Daily Treatment','FontSize',20,'FontWeight','bold')
        fname_fig = 'daily';
    elseif i == 3
        sgtitle('Adaptive Treatment','FontSize',20,'FontWeight','bold')
        fname_fig = 'adaptive';
    elseif i == 4
        sgtitle('Daily Treatment + Toxic Feedback','FontSize',20,'FontWeight','bold')
        fname_fig = 'daily_with_toxicity';
    elseif i == 5
        sgtitle('Adaptive Treatment + Toxic Feedback','FontSize',20,'FontWeight','bold')
        fname_fig = 'adaptive_with_toxicity';
    end
    
    subplot(1,2,1)
    hold on; % daily
    set(gca,'LineWidth',1.25,'FontSize',16,'FontWeight','normal','FontName','Helvetica')
    a0 = area(x,y,'edgecolor','none','FaceColor','#D3D3D3','FaceAlpha',0.5);
    h1 = plot(time_all{i},sens_all{i},'LineWidth',2,'Color',"#0072BD");
    h2 = plot(time_all{i},resist_all{i},'LineWidth',2,'Color',"#D95319");
    h3 = plot(time_all{i},sens_all{i}+resist_all{i},'LineWidth',2,'Color',"#A2142F");
    yline(Nfail,'k--','Rx Failure','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','middle')
    yline(Nhi,'k--','Rx On','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','middle')
    yline(Nlo,'k--','Rx Off','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','middle')
    legend([h1,h2,h3],' Sensitive Cells',' Resistant Cells',' Total Cells',...
        'Location','NorthEast');
    legend boxoff
    hold off;
    xlabel('Time (days)','FontSize',16)
    ylabel('Tumor Size','FontSize',16)
    xlim([0 tspan(end)])
    ylim([0 100])

    subplot(1,2,2)
    hold on;
    set(gca,'LineWidth',1.25,'FontSize',16,'FontWeight','normal','FontName','Helvetica')
    a1 = area(x,y,'edgecolor','none','FaceColor','#D3D3D3','FaceAlpha',0.5);
    g1 = plot(time_all{i},drug_all{i},'Color',"#7E2F8E",'LineWidth',2);
    g2 = plot(time_all{i},tox_all{i},'Color',"#FF0000",'LineWidth',2);
    yline(Thi,'k--','Tox Off','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','middle')
    yline(Tlo,'k--','Tox On','LineWidth',2,'FontSize',16,'LabelVerticalAlignment','middle')
    legend([g1,g2],' Drug Concentration',' Toxicity');
    legend boxoff
    hold off;
    xlabel('Time (days)','FontSize',16)
    ylabel('Drug Concentration and Toxicity','FontSize',16)
    xlim([0 tspan(end)])
    ylim([0 4])

    saveas(gcf,[fname_fig,'.fig'])
    saveas(gcf,[fname_fig,'.png'])

    clear x y
end

save output.mat tfail_all; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Updated set_parameters_ICs function
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
function [Tx, Rx, time, sens, resist, drug, tox] = solve_model(protocol, ...
    ICs, tp, nrun, q, dose, Nfail, Nlo, Nhi, Thi, Tlo)
    time = []; drug = []; tox = []; 
    sens = []; resist = []; 
    Rx = zeros(1, nrun);
    Tx = zeros(1, nrun);
    ICs(3) = dose; 
    
    fprintf('Up to protocol #%d:\n',protocol);
    on = 1;   
    for i = 1:1:nrun
        trun = (i - 1) * tp:0.1:i * tp;
        Rx(i) = on * 150;         % to draw gray background for Rx
        Tx(i) = trun(end) - 1;    % to draw gray background for Rx

        % Solve ODEs
        [t, F] = ode23s(@(t, x) tumor_model(t, x, q), trun, ICs);

        % update time and cell arrays
        time = vertcat(time, t(2:end));
        drug = vertcat(drug, F(2:end, 3));
        sens = vertcat(sens, F(2:end, 1));
        resist = vertcat(resist, F(2:end, 2));
        tox = vertcat(tox, F(2:end, 4));
    
        on_new = on; % by default, we don't change the strategy
        if protocol == 2 % daily, no toxicity feedback
            % Nothing to check, as no tumor size or toxicity feedback
        elseif protocol == 3 % adaptive, no toxicity feedback
            % check to see if tumor size condition for off cycling met
            %N = F(:, 1) + F(:, 2);
            Nend = F(end, 1) + F(end, 2);
            if on == 1 && (Nend < Nlo) %any(N < Nlo)
                on_new = 0;
                fprintf('\tTurning drug OFF at t = %f since tumor size = %f < %f\n',...
                    max(trun),Nend,Nlo);
            end
            if on == 0 && (Nend > Nhi) %any(N > Nhi)
                on_new = 1;
                fprintf('\tTurning drug ON at t = %f since tumor size = %f > %f\n',...
                    max(trun),Nend,Nhi);
            end
        elseif protocol == 4 % daily + toxicity feedback
            % check to see if toxicity condition for off cycling 
            tox_end = F(end,4); 
            if (on == 1) && (tox_end > Thi) 
                on_new = 0;
                fprintf('\tTurning drug OFF at t = %f since end(tox) = %f > %f\n',...
                    max(trun),tox_end,Thi);
            end
            if (on == 0) && (tox_end < Tlo)
                on_new = 1;
                fprintf('\tTurning drug ON at t = %f since end(tox) = %f < %f\n',...
                    max(trun),tox_end,Thi);
            end
        elseif protocol == 5 % adaptive + toxicity feedback
            Nend = F(end, 1) + F(end, 2);
            tox_end = F(end,4); 
            if (on == 1) && ((Nend < Nlo) || (tox_end > Thi))
                on_new = 0;
                fprintf('\tTurning drug OFF at t = %f for one of two reasons:\n',max(trun));
                fprintf('\t\tEither tumor size = %f < %f,\n',Nend,Nlo);
                fprintf('\t\tOr, end(tox) = %f > %f\n',tox_end,Thi);
            end
            if (on == 0) && (Nend > Nhi) && (tox_end < Tlo)
                on_new = 1;
                fprintf('\tTurning drug ON at t = %f since both:\n',max(trun));
                fprintf('\t\tTumor size = %f > %f,\n',Nend,Nhi);
                fprintf('\t\tAnd, end(tox) = %f < %f\n',tox_end,Tlo);
            end
        end
        on = on_new;

        ICs = [F(end, 1), F(end, 2), F(end, 3) + on * dose, F(end, 4)]; 
        clear t F 
        if max(sens + resist) > Nfail
            break
        end      
    end  
end