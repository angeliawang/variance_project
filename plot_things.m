% plot MM* for different weights
% single connection rate model

set(0,'DefaultAxesFontSize',15,'defaultaxeslinewidth',2,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)

% weights
weights = 0:10:100;
omega = 0:0.1:20;
expos = -1:0.01:3;
Ts = 10.^expos; %0.1:0.5:1000;
var_T = NaN(length(weights), length(Ts));
var_T_stripped = NaN(length(weights), length(Ts));
psd_T = NaN(length(weights), length(Ts));
integral_vals = NaN(length(weights), length(Ts));
M_numint = NaN(1, length(weights));
numerical_psds = NaN(length(weights), length(omega));
suffix = '_unnormed_Poiss'; % tack on at the end of filenames

% some booleans
normalize_psds = 0;
make_poiss = 0;

% other constants
zeta_0 = 0; % noise inputs to GC, zero in the model
tau = 1;
xi_0 = 1;
T = 10;
        % for some Poisson nonsense
        dt = .02;
        durationS = 1000;             % 1 sec simulation
        times = 0:dt:durationS;	% a vector with each time step	
        
        % the larger this is, the lower resolution you get at low freqs
        num_averages = 25;

fack = figure(2);
hold on
colorz = varycolor(length(weights));
for w_i = 1:length(weights)
    w = weights(w_i);
    
    syms M(x) %indep variable is omega
    M(x) = (w^2*zeta_0 + (1+x^2*tau^2)*xi_0)/((1+w)^2+x^2*(tau^2+1-2*w*tau)+x^4*tau^2);
    
    if normalize_psds
        M_int = int(M(x), x, [omega(1), omega(end)]);
        M_numint(w_i) = vpa(M_int);
        numerical_psds(w_i, :) = M(omega)/M_numint(w_i);
    else
        numerical_psds(w_i, :) = M(omega);
    end
    
    subplot(2, 1, 1)
    % this if loop nonsense is necessary becaues i only want to 
    % label the first and last 
    if w_i==1
            w_1 = plot(omega, numerical_psds(w_i, :), 'color', colorz(w_i, :));
        elseif w_i == length(weights)
            w_end = plot(omega, numerical_psds(w_i, :), 'color', colorz(w_i, :));
        else
            plot(omega, numerical_psds(w_i, :), 'color', colorz(w_i, :));
    end
    hold on
    
    for t_i = 1:length(Ts)
        T = Ts(t_i);

        % with respect to omegac
        % where
        integrand = @(x, zeta, t, xi, we, bin) (we^2*zeta + (1+x.^2*t^2)*xi)./((1+we).^2+x.^2.*(t^2+1-2*we*t)+x.^4.*t^2).*(1-cos(x*bin))./x.^2;
%         psd_bit = @(x, zeta, t, xi, we, bin) (we^2*zeta + (1+x.^2*t^2)*xi)./((1+we).^2+x.^2.*(t^2+1-2*we*t)+x.^4.*t^2);
        var_T(w_i, t_i) = integral(@(x) integrand(x,zeta_0,tau,xi_0,w,T), -Inf, Inf)/((pi*T^2));
        
        if normalize_psds
            subplot(2, 1, 2)
            if w_i==1
                w2_1 = loglog(Ts, var_T(w_i, :)./M_numint(w_i), 'color', colorz(w_i, :));
            elseif w_i==length(weights)
                w2_end = loglog(Ts, var_T(w_i, :)./M_numint(w_i), 'color', colorz(w_i, :));
            else
                loglog(Ts, var_T(w_i, :)./M_numint(w_i), 'color', colorz(w_i, :))
            end
            hold on
        else
            subplot(2, 1, 2)
            if w_i==1
                w2_1 = loglog(Ts, var_T(w_i, :), 'color', colorz(w_i, :));
            elseif w_i==length(weights)
                w2_end = loglog(Ts, var_T(w_i, :), 'color', colorz(w_i, :));
            else
                loglog(Ts, var_T(w_i, :), 'color', colorz(w_i, :))
            end
            hold on
        end
    end
    
    if make_poiss
        % first we have the variances of rates, mult by the corresponding
        % window size to get the spike counts 
        mean_spikecounts = var_T(w_i, :).*Ts;
        
        Ts_floor = NaN(size(Ts));
        varsim_T = NaN(1, length(Ts));
        
        for tsize = 1:length(Ts)
             % now M_numint is the mean variance.
            exp_rate = mean_spikecounts; % rate in HZ?
            spike_probs = rand(size(times));
            poiss_spikes = (exp_rate*dt) > spike_probs;
        
            if w_i==1
                [omega, psds] = power_spectrum(poiss_spikes, dt, num_averages);
                numerical_psds_poiss = NaN(length(weights), length(omega));
            end
            numerical_psds_poiss(w_i, :) = power_spectrum(poiss_spikes, dt, num_averages); 
            
        % now compute the variance by taking a buncha samples of varying T
        num_samples = 200;
        vsim_fig = figure;

        for t_i = 1:length(Ts)
            Tf = ceil(Ts(t_i));
            Ts_floor(t_i) = Tf;
            sample_length = Tf/dt-1;
            num_steps = Tf/dt;
            samples = NaN(num_samples, sample_length+1);

            % populate the samples matrix
            for ns_i = 1:num_samples
                start = randi(num_steps-sample_length);
                samples(ns_i, :) = poiss_spikes(start:start+sample_length);
            end
            % this is spiking data

            % take the sum and div by T and take the var      
            varsim_T(t_i) = nanvar(sum(samples, 2)/Tf);
        end

        subplot(2, 1, 2)
        if w_i==1
        w2_1 = loglog(Ts_floor, varsim_T, 'color', colorz(w_i, :));
        elseif w_i==length(weights)
        w2_end = loglog(Ts_floor, varsim_T, 'color', colorz(w_i, :));
        else
            loglog(Ts_floor, varsim_T, 'color', colorz(w_i, :));
        end
        hold on 
            
        end
        
    end
end
hold off

subplot(2, 1, 1)
set(gca, 'Yscale', 'log')
xlabel('\omega')
ylabel('\mid M(\omega) \mid^2')
legend([w_1, w_end], ['w=', num2str(weights(1))], ...
    ['w=', num2str(weights(end))], 'location', 'northeast')
subplot(2, 1, 2)
xlabel('T', 'interpreter', 'latex')
ylabel('var ( $\frac{1}{T} \int_0^T M(t) dt$ )', 'interpreter', 'latex')
legend([w2_1, w2_end], ['w=', num2str(weights(1))], ...
    ['w=', num2str(weights(end))], 'location', 'southwest')
% subplot(3, 1, 3)
% xlabel('T', 'interpreter', 'latex')
% ylabel('Without normalization by T^2', 'interpreter', 'latex')
% legend([w3_1, w3_end], 'w=0', 'w=100', 'location', 'southwest')
hgsave(fack, ['psd_and_var', suffix,'.fig'])
saveas(fack, ['psd_and_var', suffix,'.png'])

% now we save some shite:
save(['var_psd_data', suffix,'.mat'], 'var_T', 'weights', ...
    'Ts', 'omega', 'M_numint', 'numerical_psds')

% %% now we look at why there are bumps at certain values of T
% diagnostic_indices = [1, 3, 6];
% % diagnostic_var = var_T(1, :);
% % diagnostic_var2 = var_T(3, :);
% diagnostic_var3 = var_T(6, :);
% 
% % % includes the cosine term how
% integrand = @(x, zeta, t, xi, we, bin) (we^2*zeta + (1+x.^2*t^2)*xi)./...
%             ((1+we).^2+x.^2.*(t^2+1-2*we*t)+x.^4.*t^2).*(1-cos(x*bin))./x.^2;
% 
% % does not include the cosine term
% psd_bit = @(x, zeta, t, xi, we, bin) (we^2*zeta + (1+x.^2*t^2)*xi)./...
%     ((1+we).^2+x.^2.*(t^2+1-2*we*t)+x.^4.*t^2);
% 
% % % gets integrated over omega 
% % % perhaps we'll just pick one and look at T
% % 
% % % or make a heatmap across different ones
% oms = -16.1:0.2:16.1;
% % % Ts is already defined
% var_map = NaN(length(oms), length(Ts));
% % psd_map = NaN(length(oms), length(Ts));
% % var_map2 = NaN(length(oms), length(Ts));
% % psd_map2 = NaN(length(oms), length(Ts));
% var_map3 = NaN(length(oms), length(Ts));
% % psd_map3 = NaN(length(oms), length(Ts));
% % 
% for om_i = 1:length(oms)
%     for t_i = 1:length(Ts)
%         % equals 1 at omega = 0 as expected
%         var_map(om_i, t_i) = integrand(oms(om_i), 1, 1, 1, 0, Ts(t_i));
% %         psd_map(om_i, t_i) = psd_bit(oms(om_i), 1, 1, 1, 0, Ts(t_i));
%         
%         % nonzero weights
%         % now there is a divot near omega = 0
% %         var_map2(om_i, t_i) = integrand(oms(om_i), 1, 1, 1, 2, Ts(t_i));
% %         psd_map2(om_i, t_i) = psd_bit(oms(om_i), 1, 1, 1, 2, Ts(t_i));
%         var_map3(om_i, t_i) = integrand(oms(om_i), 1, 1, 1, 5, Ts(t_i));
% %         psd_map3(om_i, t_i) = psd_bit(oms(om_i), 1, 1, 1, 5, Ts(t_i));
%     end
% end
% 
% sum_across_oms = squeeze(sum(var_map3, 1));
% sums_across_oms = repmat(sum_across_oms, length(oms), 1);
% normed_var_map3 = var_map3./sums_across_oms;
% 
% figure; 
% imagesc(Ts, oms, normed_var_map3)
% ylabel('\omega')
% xlabel('T')
% 
% % title('w = 0')
% % subplot(3, 1, 2)
% % imagesc(Ts, oms, var_map2)
% % title('w = 2')
% % subplot(3, 1, 3)
% % imagesc(Ts, oms, var_map3)
% % title('w = 5')
% % 
% % figure;
% % for d_i = 1:length(diagnostic_indices)
% % loglog(Ts, var_T(d_i, :)/M_numint(d_i))
% % hold on
% % end
% % hold off
% % xlabel('T')
% % ylabel('var')
% % 
% % % maybe it's worth looking at the T values
% 
% dist = @(x, bin) (1-cos(x*bin))./x.^2./(pi*bin);
% 
% omega_expos = -4:0.001:2;
% omegas = [0, 10.^omega_expos];
% trial_expos = -1:0.5:2;
% trial_Ts = 10.^trial_expos;
% T_colorz = varycolor(length(trial_Ts));
% 
% cos_fig = figure;
% for tT_i = 1:length(trial_Ts)
%     trial_T = trial_Ts(tT_i);
%     if tT_i==1 
%         cos_1 = plot(omegas, dist(omegas, trial_T), 'color', T_colorz(tT_i,:));
%         hold on
%     elseif tT_i==length(trial_Ts)
%         cos_end = plot(omegas, dist(omegas, trial_T), 'color', T_colorz(tT_i,:));
%         hold on
%     else
%         plot(omegas, dist(omegas, trial_T), 'color', T_colorz(tT_i,:));
%         hold on
%     end
% end
% set(gca, 'Yscale', 'log')
% set(gca, 'Xscale', 'log')
% xlabel('$\omega$', 'interpreter', 'latex')
% legend([cos_1, cos_end], 'T=0.1', 'T=1000', 'location', 'northeast')
% hgsave(cos_fig, 'cos_bit.fig')
% saveas(cos_fig, 'cos_bit.png')