%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                       %
% Quantifying the impact of treatment-related thresholds across a       %
% virtual population.                                                   %
% Updated: 2/24/2025                                                    %
%                                                                       %
% This code expands upon sweep_VP_5params by varying ALL FOUR           %
% thresholds simultaneously to assesses VP response the adaptive +      %
% toxicity protocol. This is the only protocol in which all four        %
% thresholds influence treatment response.                              %
%                                                                       %
% Run time desktop for 100 VPs: 14496 sec = 241 min.6 min ~= 4 hrs      %
%                                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc; tic;
rng(1); % fixed seed for replicability
[q, ICs] = set_parameters_ICs();
num_VPs = 100; % Number of virtual patients 
num_pts = 13; % Discretization of protocol thresholds 

%% VP distributions: Parameter means and standard deviations
% param 1 = alphaS (range 0.5-1.5), 2 = epsilon (0.2-1), 
% 3 = beta (0.5-3.5), 4 = delta (0.5-1.5), 5 = gamma (0.2-0.8)
param_mean = [q.alphaS q.epsilon q.beta q.delta q.gamma]; 
param_std  = [0.18     0.12      0.45   0.16    0.09];
param_labels = {'\alpha_S', '\epsilon', '\beta', '\delta', '\gamma'};
VP_params = generate_vp_params(num_VPs, param_mean, param_std);
min_val = min(VP_params);
max_val = max(VP_params); 
range_val = max_val-min_val;
x_limits = zeros(length(param_labels), 2);
x_limits(:,1) = min_val' - 0.3 * range_val';
x_limits(:,2) = max_val' + 0.3 * range_val';

figure;
set(gcf, 'Units', 'Normalized','OuterPosition', [0.05, 0.05, 0.5, 0.7]);
sgtitle('Distributions of Virtual Patient Parameters', 'FontSize', 18, 'FontWeight', 'bold');
for i = 1:length(param_labels)
    % Plot VP parameters
    subplot(2,3, i);
    histogram(VP_params(:, i), 'Normalization', 'pdf');
    hold on;
    
    % Plot the theoretical lognormal PDF for each parameter
    mu_Y = log(param_mean(i)) - ...
        0.5 * log((param_std(i)^2 / param_mean(i)^2) + 1);
    sigma_Y = sqrt(log((param_std(i)^2 / param_mean(i)^2) + 1));
    
    % Generate x-values for the PDF plot
    x = linspace(x_limits(i, 1), x_limits(i, 2), 1000);
    pdf_theoretical = lognpdf(x, mu_Y, sigma_Y);
    plot(x, pdf_theoretical, 'r', 'LineWidth', 2);
    
    xlabel(param_labels{i}, 'FontSize', 16);
    %ylabel('Probability Density', 'FontSize', 16);
    title(['\mu = ' num2str(param_mean(i)) ', \sigma = ' num2str(param_std(i))], ...
        'FontSize', 18);
    xlim(x_limits(i, :));

    % Display range of each parameter
    fprintf('%s: ',param_labels{i});
    disp([min(VP_params(:, i)) max(VP_params(:, i))])
end

%% Protocol settings
tspan = 0:0.1:100;
nrun = max(tspan);           % number of therapy cycles for continuous
tp = 1;                      % therapy period
dose = 1;                    % drug bolus dose
baseline = (ICs(1) + ICs(2));
Nhi = 1.0*baseline;   % adaptive threshold: Rx turns on
Nlo = 0.4*baseline;   % adaptive threshold: Rx turns off 
Nfail = 1.5 * baseline; % Rx failure threshold
Thi = 2;                     % Toxicity tolerance: upper bound
Tlo = 1;                     % Toxicity tolerance: lower bound
Tlo_vec = linspace(0.25, 2.25, num_pts); 
Thi_vec = linspace(0.5, 4, num_pts); 
Nhi_vec = linspace(25, Nfail, num_pts); % must keep above Nlo = 20 and below Nfail
Nlo_vec = linspace(5, 45, num_pts); 
protocol_string = {'', 'Daily', 'Adaptive', 'Daily + Toxicity', 'Adaptive + Toxicity'};

%% Get average behavior of VPs in scenarios where thresholds do matter
alphaS_VPs = VP_params(:,1);
epsilon_VPs = VP_params(:,2);
beta_VPs = VP_params(:,3);
delta_VPs = VP_params(:,4);
gamma_VPs = VP_params(:,5);

%%  Vary all 4 protocol thresholds for adaptive + toxicity
protocol = 5; 
agg_tfail_sweep_all = zeros(length(Tlo_vec), length(Thi_vec),length(Nhi_vec),...
    length(Nlo_vec),num_VPs); 
for vp = 1:num_VPs
    % Set the parameters for this virtual patient
    q.beta = beta_VPs(vp);
    q.gamma = gamma_VPs(vp);
    q.alphaS = alphaS_VPs(vp);
    q.delta = delta_VPs(vp);
    q.epsilon = epsilon_VPs(vp);
    q.alphaR = q.epsilon * q.alphaS; % Using epsilon

    % Display progress and individual values
    fprintf('Running 4D protocol sweep for VP #%d of %d:\n', vp,num_VPs);
    fprintf('\tbeta: %.2f, gamma: %.2f, alphaS: %.2f, delta: %.2f, epsilon: %.2f\n', ...
        q.beta, q.gamma, q.alphaS, q.delta, q.epsilon);

     %% Run the analysis for this patient for only adaptive + tox (protocol 5)
    tfail = zeros(length(Tlo_vec), length(Thi_vec),length(Nhi_vec), length(Nlo_vec)); 
    for i = 1:length(Tlo_vec)
        fprintf('\t\tUp to Tlo = %f\n',Tlo_vec(i));
        for j = 1:length(Thi_vec)
           for k = 1:length(Nhi_vec)
                for l = 1:length(Nlo_vec)
                    if Thi_vec(j) <= Tlo_vec(i)
                        tfail(i,j,k,l) = nan;
                    elseif Nhi_vec(k) <= Nlo_vec(l)
                        tfail(i,j,k,l) = nan;
                    else % valid protocol
                        [Tx, Rx, time, sens, resist, drug, tox] = ...
                            solve_model(protocol, ICs, tp, nrun, q, dose, ...
                            Nfail, Nlo_vec(l), Nhi_vec(k), Thi_vec(j), Tlo_vec(i));
                        time_to_check = time(find((sens + resist) > Nfail, 1));
                        if isempty(time_to_check)
                            time_to_check = max(tspan) * 1.5;
                        end
                        tfail(i,j,k,l) = time_to_check;
                    end
                end
           end
        end
    end
    agg_tfail_sweep_all(:,:,:,:,vp) = tfail;
end
toc; 

% Calculate the mean values across VPs in scenarios where thresholds do matter
mean_tfail_sweep_all = nanmean(agg_tfail_sweep_all, 5);
max_val_all = max(mean_tfail_sweep_all,[],'all');

% Extract index of optimal protocol
maxIdx_lin = find(mean_tfail_sweep_all==max_val_all);
[Tlo_maxIdx, Thi_maxIdx, Nhi_maxIdx, Nlo_maxIdx] = ...
    ind2sub(size(mean_tfail_sweep_all),maxIdx_lin);

fprintf('Adaptive+tox has max avg TTP of %f with:\n',max_val_all);
fprintf('\tToxicity thresholds: Tlo = %f, Thi = %f\n',...
    Tlo_vec(Tlo_maxIdx),Thi_vec(Thi_maxIdx));
fprintf('\tAdaptive thresholds: Relative Nhi = %f, relative Nlo = %f\n',...
    Nhi_vec(Nhi_maxIdx)/baseline,Nlo_vec(Nlo_maxIdx)/baseline);

save Output_VPs/output_VPs_4D.mat param_mean param_std VP_params Tlo_vec ...
    Thi_vec Nhi_vec Nlo_vec protocol_string agg_tfail_sweep_all ....
    mean_tfail_sweep_all max_val_all Tlo_maxIdx Thi_maxIdx Nhi_maxIdx Nlo_maxIdx

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
function [Tx, Rx, time, sens, resist, drug, tox] = solve_model(protocol, ...
    ICs, tp, nrun, q, dose, Nfail, Nlo, Nhi, Thi, Tlo)
    time = []; drug = []; tox = []; 
    sens = []; resist = []; 
    Rx = zeros(1, nrun);
    Tx = zeros(1, nrun);
    ICs(3) = dose; 
  
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
      
        ICs = [F(end, 1), F(end, 2), F(end, 3) + on * dose, F(end, 4)]; 
        clear t F 
        if max(sens + resist) > Nfail
            break
        end      
    end  
end

function VP_params = generate_vp_params(N_VPs, param_mean, param_std)
    VP_params = zeros(N_VPs, length(param_mean));
    for i = 1:length(param_mean)
        mu_X = param_mean(i);
        sigma_X = param_std(i);
        
        sigma_Y_sq = log((sigma_X^2 / mu_X^2) + 1);
        sigma_Y = sqrt(sigma_Y_sq);
        mu_Y = log(mu_X) - (sigma_Y_sq / 2);
        
        VP_params(:, i) = lognrnd(mu_Y, sigma_Y, [N_VPs, 1]);
    end
end
