% plot MM* for different weights
% single connection rate model

set(0,'DefaultAxesFontSize',15,'defaultaxeslinewidth',2,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)

% weights
weights = 0:1:10;
omega = 0:0.1:20;
expos = -3:0.05:3;
Ts = 10.^expos; %0.1:0.5:1000;

covar_T = NaN(length(weights), length(Ts));
psd1_T = NaN(length(weights), length(Ts));
psd2_T = NaN(length(weights), length(Ts));

var_T_stripped = NaN(length(weights), length(Ts));
integral_vals = NaN(length(weights), length(Ts));
M_numint = NaN(1, length(weights));
numerical_psds = NaN(length(weights), length(omega));

% which variable are we changing?
indep_var = weights;
num_vars = length(indep_var);
indep_xlabel = 'weights';

suffix = ['_2MC_z1']; % tack on at the end of filenames

% some booleans
normalize_xcorr = 0; % if we want the integral to all be 1
make_poiss = 0;

% other constants
zeta_0 = 1; % noise inputs to GC, zero in the model
tau = 1;
xi_1_0 = 1;
xi_2_0 = 1;

if xi_1_0==xi_2_0
    xi_0 = xi_1_0;
    fprintf('Noise distributions are the same \n')
else
    fprintf('Noise distributions are not the same \n')
end
T = 10;

colorz = varycolor(num_vars);

integrand = @(x, zeta, t, xi1, xi2, we, bin) ...
    ((-we.*xi1.*(1-x.^2*t+we)-we.*xi2.*(1-x.^2.*t+we)+we.^2.*zeta.*(1+x.^2))./...
((1+x.^2).*((1+x.^2).*(1+x.^2.*t^2)+4.*we.*(1-x.^2.*t)+4.*we.^2))).*(1-cos(x.*bin))./x.^2;

    % diagonal ones
integrand2 = @(x, zeta, t, xi, we, bin) ...
    (xi*((1+x.^2).*(1+x.^2.*t^2)+2*we*(1-x.^2.*t)+2*we^2) + we^2.*zeta.*(1+x.^2))./...
((1+x.^2).*((1+x.^2).*(1+x.^2.*t^2)+4.*we.*(1-x.^2.*t)+4.*we.^2)).*(1-cos(x.*bin))./x.^2;

fack = figure(2);
hold on
for w_i = 1:num_vars
    %% MAY NEED TO BE ADJUSTED, depending on what the var is
    w = indep_var(w_i);
    
    % this is the OFF-DIAGONAL CROSS-CORR term.
    syms M(x) %indep variable is omega
    M(x) = (-w*xi_1_0*(1-x.^2*tau+w)-w*xi_2_0*(1-x.^2*tau+w)+w.^2*zeta_0*(1+x.^2))/...
        ((1+x.^2)*((1+x.^2)*(1+x.^2*tau^2)+4*w*(1-x.^2*tau)+4*w.^2));
    
    % these are the DIAGONAL PSD terms
    syms P1(x)
    P1(x) = (xi_1_0*((1+x.^2)*(1+x.^2.*tau^2)+2*w*(1-x.^2*tau)+w^2) + w^2*xi_2_0 + w^2*zeta_0*(1+x.^2))./...
        ((1+x.^2)*((1+x.^2)*(1+x.^2*tau^2)+4*w*(1-x.^2*tau)+4*w.^2));
    
    syms P2(x)
    P2(x) = (xi_2_0*((1+x.^2)*(1+x.^2.*tau^2)+2*w*(1-x.^2*tau)+w^2) + w^2*xi_1_0 + w^2*zeta_0*(1+x.^2))./...
        ((1+x.^2)*((1+x.^2)*(1+x.^2*tau^2)+4*w*(1-x.^2*tau)+4*w.^2));
    
    if normalize_xcorr
        % find the integral
        M_int = int(M(x), x, [omega(1), omega(end)]);
        M_numint(w_i) = vpa(M_int);
        numerical_psds(w_i, :) = M(omega)/M_numint(w_i);
    else
        numerical_psds(w_i, :) = M(omega);
    end
    
    subplot(2, 3, 1)
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

        % with respect to omega
        % where x = omega
        % t = tau
        % bin = T
    
        zeta = zeta_0; t = tau; xi1 = xi_1_0; xi2 = xi_2_0; we = w; bin = T;
        
        covar_T(w_i, t_i) = integral(@(x) ((-we.*xi1.*(1-x.^2*t+we)-we.*xi2.*(1-x.^2.*t+we)+we.^2.*zeta.*(1+x.^2))./...
        ((1+x.^2).*((1+x.^2).*(1+x.^2.*t^2)+4.*we.*(1-x.^2.*t)+4.*we.^2))).*(1-cos(x.*bin))./x.^2, -Inf, Inf)/(pi*bin^2);
        
        psd1_T(w_i, t_i) = integral(@(x) (xi1*((1+x.^2).*(1+x.^2.*t^2)+2*we*(1-x.^2.*t)+we^2) + we.^2*xi2 + we^2.*zeta.*(1+x.^2))./...
        ((1+x.^2).*((1+x.^2).*(1+x.^2.*t^2)+4.*we.*(1-x.^2.*t)+4.*we.^2)).*(1-cos(x.*bin))./x.^2, -Inf, Inf)/(pi*bin^2);
        
        psd2_T(w_i, t_i) = integral(@(x) (xi2*((1+x.^2).*(1+x.^2.*t^2)+2*we*(1-x.^2.*t)+we^2) + we.^2*xi1 + we^2.*zeta.*(1+x.^2))./...
        ((1+x.^2).*((1+x.^2).*(1+x.^2.*t^2)+4.*we.*(1-x.^2.*t)+4.*we.^2)).*(1-cos(x.*bin))./x.^2, -Inf, Inf)/(pi*bin^2);
    end
%     integrand2(3,zeta_0,tau,xi_0,w,T)-integrand(3,zeta_0,tau,xi_1_0,xi_2_0,w,T)
    
    if normalize_xcorr
        subplot(2, 3, 4)
        if w_i==1
            w2_1 = loglog(Ts, covar_T(w_i, :)./M_numint(w_i), 'color', colorz(w_i, :));
        elseif w_i==length(weights)
            w2_end = loglog(Ts, covar_T(w_i, :)./M_numint(w_i), 'color', colorz(w_i, :));
        else
            loglog(Ts, covar_T(w_i, :)./M_numint(w_i), 'color', colorz(w_i, :))
        end
        hold on
    else
        subplot(2, 3, 4)
        if w_i==1
            w2_1 = loglog(Ts, covar_T(w_i, :), 'color', colorz(w_i, :));
        elseif w_i==length(weights)
            w2_end = loglog(Ts, covar_T(w_i, :), 'color', colorz(w_i, :));
        else
            loglog(Ts, covar_T(w_i, :), 'color', colorz(w_i, :))
        end
        hold on
    end

    subplot(2, 3, 2)
    if w_i==1
        w3_1 = loglog(Ts, psd1_T(w_i, :), 'color', colorz(w_i, :));
    elseif w_i==length(weights)
        w3_end = loglog(Ts, psd1_T(w_i, :), 'color', colorz(w_i, :));
    else
        loglog(Ts, psd1_T(w_i, :), 'color', colorz(w_i, :))
    end
    hold on
    
    subplot(2, 3, 5)
    if w_i==1
        w6_1 = loglog(Ts, psd2_T(w_i, :), 'color', colorz(w_i, :));
    elseif w_i==length(weights)
        w6_end = loglog(Ts, psd2_T(w_i, :), 'color', colorz(w_i, :));
    else
        loglog(Ts, psd2_T(w_i, :), 'color', colorz(w_i, :))
    end
    hold on

    subplot(2, 3, 3)
    if w_i==1
        w4_1 = loglog(Ts, psd1_T(w_i, :)-covar_T(w_i, :), 'color', colorz(w_i, :));
    elseif w_i==length(weights)
        w4_end = loglog(Ts, psd1_T(w_i, :)-covar_T(w_i, :), 'color', colorz(w_i, :));
    else
        loglog(Ts, psd1_T(w_i, :)-covar_T(w_i, :), 'color', colorz(w_i, :))
    end
    hold on
    
    subplot(2, 3, 6)
    if w_i==1
        w5_1 = loglog(Ts, psd1_T(w_i, :)+covar_T(w_i, :), 'color', colorz(w_i, :));
    elseif w_i==length(weights)
        w5_end = loglog(Ts, psd1_T(w_i, :)+covar_T(w_i, :), 'color', colorz(w_i, :));
    else
        loglog(Ts, psd1_T(w_i, :)+covar_T(w_i, :), 'color', colorz(w_i, :))
    end
    hold on
    
end
hold off

figure(2);
subplot(2, 3, 1)
set(gca, 'yscale', 'linear')
xlabel('\omega')
ylabel('$\frac{1}{2} (M_1^* M_2 + M_2 M_1^*)$', 'interpreter', 'latex')
axis square;
legend([w_1, w_end], ['w=', num2str(weights(1))], ...
    ['w=', num2str(weights(end))], 'location', 'southeast')

subplot(2, 3, 4)
set(gca, 'yscale', 'linear')
xlabel('T', 'interpreter', 'latex')
ylabel('$Cov(M_1, M_2)$', 'interpreter', 'latex')
%ylabel('var ( $\frac{1}{T} \int_0^T M_1(t_1) dt_1 \frac{1}{T} \int_0^T M_2(t_2) dt_2 $ )', 'interpreter', 'latex')
axis square;
legend([w2_1, w2_end], ['w=', num2str(weights(1))], ...
    ['w=', num2str(weights(end))], 'location', 'southeast')

subplot(2, 3, 2)
set(gca, 'yscale', 'linear')
xlabel('T', 'interpreter', 'latex')
ylabel('$ Var(M_1)) $', 'interpreter', 'latex')
axis square;
legend([w3_1, w3_end], ['w=', num2str(weights(1))], ...
    ['w=', num2str(weights(end))], 'location', 'northeast')

subplot(2, 3, 5)
set(gca, 'yscale', 'linear')
xlabel('T', 'interpreter', 'latex')
ylabel('$ Var(M_2)$', 'interpreter', 'latex')
axis square;
legend([w6_1, w6_end], ['w=', num2str(weights(1))], ...
    ['w=', num2str(weights(end))], 'location', 'northeast')

subplot(2, 3, 3)
set(gca, 'yscale', 'linear')
xlabel('T', 'interpreter', 'latex')
ylabel('odd $\lambda$', 'interpreter', 'latex')
axis square;
legend([w4_1, w4_end], ['w=', num2str(weights(1))], ...
    ['w=', num2str(weights(end))], 'location', 'northeast')

subplot(2, 3, 6)
set(gca, 'yscale', 'linear')
xlabel('T', 'interpreter', 'latex')
ylabel('even $\lambda$', 'interpreter', 'latex')
axis square;
legend([w5_1, w5_end], ['w=', num2str(weights(1))], ...
    ['w=', num2str(weights(end))], 'location', 'northeast')

hgsave(fack, ['variance_studies/psd_var', suffix,'.fig'])
saveas(fack, ['variance_studies/psd_var', suffix,'.png'])

% now we save some shite:
if normalize_xcorr
    save(['variance_studies/var_psd_data', suffix,'.mat'], 'var_T', 'weights', ...
    'Ts', 'omega', 'M_numint', 'numerical_psds')
else
    save(['variance_studies/var_psd_data', suffix,'.mat'], 'covar_T', 'weights', 'Ts', 'omega')
end
