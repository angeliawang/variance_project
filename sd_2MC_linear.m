dirs = {'2MC_single_linear_fastGC/', '2MC_global_linear_fastGC/'};
suffix = '';

set(0, 'DefaultFigureVisible', 'on')
set(0,'DefaultAxesFontSize',30,'defaultaxeslinewidth',2,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)

% if you've already computed all the information and just want to plot stuff
just_plotting = 0;

num_copies = 20;
inhib_strengths = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 5, 10, 20, 50]; % for all types of connectivity
windows = 5:5:200;
Nm = 2;
var_of_interest = inhib_strengths;
xlabel_of_interest = 'inhibition';

am_i_on_quest = 1;

dt = .02;
T_final = 20000;
Nt = T_final/dt;

% number of data points for LDA
% i.e. number of random windows to take
num_samples = 200;

runs = 20;

assert(num_samples>= Nm);
% else the covariance matrix is guaranteed singular

% number of times to run this, to smooth things out

for reg = dirs
        
    target_dir = ['AdExIF/', reg{1}];
    files = dir(target_dir);
    directoryNames = {files([files.isdir]).name};
    directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

    if ~just_plotting
    % usually these are the same
    % there is the odd run where quest timed out and we didn't hit all the
    % inhib weights
    sup = min(length(var_of_interest), length(directoryNames)/num_copies);
    
    % initialize things
    % we use zeros() instead of NaN() because I use a rolling average
    F_opt = zeros(sup, length(windows), num_copies, runs); % one for Poiss and one for netw
    pc_opt = zeros(size(F_opt));
    mc_frs_1 = zeros(Nm, sup, num_copies);
    mc_frs_2 = zeros(Nm, sup, num_copies);
    
    covsums = NaN(Nm, Nm, sup, length(windows), num_copies, runs);
    covsum_evalues = zeros(Nm, sup, length(windows), num_copies, runs);
    covsum_evecs = NaN(Nm, Nm, sup, length(windows), num_copies, runs);

        corrcoefs = zeros(Nm, Nm, sup, length(windows), num_copies, runs);
        diff_mean_vector = zeros(Nm, sup, length(windows), num_copies, runs);
        Fopt_mode = NaN(size(covsum_evalues));
        Fopt_mc = NaN(size(Fopt_mode));
        
    for d_i = 1:size(directoryNames, 2) 
	try

        d_i
        file_location = [target_dir,directoryNames{d_i}];

        full_path = [file_location, '/data.mat'];
    %     full_path = [file_location, '/', file.name];
        data_mc = load(full_path);
        
        Voltage = data_mc.Voltage_history;
        %Nm = size(Voltage, 1);
        Nm = 2;
       Voltage = Voltage(1:Nm, :, :);

            
        Ns = size(Voltage, 2);
        num_steps = size(Voltage, 3);

        % Do some stuff

        % inhib_strength indices
        i_i = ceil(d_i/num_copies);
        ii_i = mod(d_i, num_copies); % from 1 to num_copies
        if ii_i==0
            ii_i = num_copies;
        end

        % mc_frs are Nm by num_inhibs by num_copies
        % and contain *FIRING RATES*
        mc_frs_1(:, i_i, ii_i) = sum(Voltage(1:Nm, 1, :), 3)/(T_final); 
        mc_frs_2(:, i_i, ii_i) = sum(Voltage(1:Nm, 2, :), 3)/(T_final); 

        try 
        %% for each wdow
        for w_i = 1:length(windows)
            wdow = windows(w_i);
            window_steps = round(wdow/dt)-1;
            
            % set the seed for debugging
%             rng(w_i);

                data_set = NaN(Nm, Ns, num_samples);
                
                for r_i = 1:runs

                % for each sample, grab that chunk
                for i = 1:num_samples

                    i1 = randi(num_steps-window_steps);
                    i2 = randi(num_steps-window_steps);

                    % these contain SPIKE COUNTS in a designated wdow
                    data_set(:, 1, i) = sum(Voltage(:, 1, i1:i1+window_steps), 3);
                    data_set(:, 2, i) = sum(Voltage(:, 2, i2:i2+window_steps), 3);
                end

%                fprintf('we are here\n')
                stim1 = squeeze(data_set(:, 1, :));
                stim2 = squeeze(data_set(:, 2, :));

                cc = corrcoef(stim1') + corrcoef(stim2');
                % corrcoefs is Nm x Nm x inhibs x windows x
                % (num_copies*num_runs)
                corrcoefs(:, :, i_i, w_i, (ii_i-1)*runs + r_i) = cc;


%                 if ii_i==1
%                     sample_data(2, 1, 1:length(channel_keep_mc), :, i_i, w_i) = ...
%                         stim1;
%                     sample_data(2, 2, 1:length(channel_keep_mc), :, i_i, w_i) = ...
%                         stim2;
%                 end
                
                % these are Nm by 1
                stim1_mean = mean(stim1, 2)/wdow;
                stim2_mean = mean(stim2, 2)/wdow;

                % in kHz
                stim1 = stim1/wdow;
                stim2 = stim2/wdow;

                % netw
                diff_means = stim1_mean-stim2_mean;
                cov_sum = cov(stim1')+cov(stim2');
                w = cov_sum\diff_means;
                c = dot(w, (stim1_mean+stim2_mean))/2;
                P1c = w'*stim1;
                P2c = w'*stim2;
                pc_opt(i_i, w_i, (ii_i-1)*runs + r_i) = ...
                    sum([P1c>c P2c<c])/num_samples/2;
                F_opt(i_i, w_i, (ii_i-1)*runs + r_i) = (dot(w, diff_means))^2/(w'*cov_sum*w);
                
                covsums(:, :, i_i, w_i, (ii_i-1)*runs + r_i) = cov_sum;

                [vecs, vals] = eig(cov_sum);
                [vals_s, ind] = sort(diag(vals));
                vecs = vecs(:, ind);
                covsum_evalues(:, i_i, w_i, ...
                    (ii_i-1)*runs + r_i) = vals_s;
                covsum_evecs(:, :, i_i, w_i, (ii_i-1)*runs + r_i) = vecs;
                
                for mc_i = 1:Nm
                    Fopt_mode(mc_i, i_i, w_i, (ii_i-1)*runs + r_i) = ...
                        dot(vecs(:, mc_i), diff_means)^2/vals_s(mc_i);
                    Fopt_mc(mc_i, i_i, w_i, (ii_i-1)*runs + r_i) = ...
                        diff_means(mc_i)*w(mc_i);
                end

                diff_mean_vector(:, i_i, w_i, (ii_i-1)*runs + r_i) = diff_means;
            end
        end
        catch
        continue
        end

catch
continue
end

end

    % corrcoefs is Nm x Nm x inhibs x windows x num_copies
    % covsum_evalues is Nm x inhibs x windows x num_copies
    
%     if runs > 1
%     F_opt = squeeze(nanmean(F_opt, 5));
%     pc_opt = squeeze(nanmean(F_opt, 5));
%     end
    
    save(strcat(target_dir, 'LDA.mat'), 'F_opt', 'pc_opt', 'mc_frs_1', ...
        'mc_frs_2')
%         save(strcat(target_dir, 'Spectra.mat'), 'covsum_evalues', 'covsum_evecs')
        save(strcat(target_dir, 'Frills.mat'), 'corrcoefs', 'diff_mean_vector', ...
            'covsum_evalues', 'covsum_evecs', 'covsums', 'Fopt_mode', 'Fopt_mc')
        
    else
        %% we already have made the things and now we just plottin
        LDA = load([target_dir, 'LDA', suffix, '.mat']);
        mc_frs_1 = LDA.mc_frs_1;
        mc_frs_2 = LDA.mc_frs_2;
        F_opt = LDA.F_opt;
        pc_opt = LDA.pc_opt;
        
        frills = load([target_dir, 'Frills', suffix, '.mat']);
        diff_mean_vector = frills.diff_mean_vector;
        mean_diffmean = squeeze(nanmean(squeeze(nanmean(diff_mean_vector, 5)), 4));
        covsums = frills.covsums;
        mean_covsums = squeeze(nanmean(squeeze(nanmean(covsums, 6)), 5)); % Nm x Nm x inhib x windows
        
        %% plotting 
        % getcha plotting here
        window_colorz = varycolor(length(windows));
        inhib_colorz = varycolor(length(inhib_strengths));
        
        % inhib strengths
            fr_1 = nanmean(mc_frs_1, 3);
            fr_2 = nanmean(mc_frs_2, 3);
            std_1 = nanstd(mc_frs_1, 0, 3);
            std_2 = nanstd(mc_frs_2, 0, 3);
        
        fr_fig = figure;
        hold on
        for i = 1:length(inhib_strengths)
            plot(fr_1(:, i), 'color', inhib_colorz(i, :))
            plot(fr_2(:, i), '--', 'color', inhib_colorz(i, :))
        end
        ylabel('average value')
        xlim([0.5, 2.5])
        hold off
        hgsave(fr_fig, [target_dir, 'mc_frs.fig'])
        saveas(fr_fig, [target_dir, 'mc_frs.png'])
    
%         fr_fig2 = figure;
%         hold on
%         for i = 1:length(inhib_strengths)
%             errorbar(fr_1(:, i), std_1(:, i), 'color', colorz(i, :))
%             errorbar(fr_2(:, i), std_2(:, i), '--', 'color', colorz(i, :))
%         end
%         ylabel('kHz')
%         xlim([0.5, 2.5])
%         hold off
%         hgsave(fr_fig2, [target_dir, 'mc_frs_errorbars.fig'])
%         saveas(fr_fig2, [target_dir, 'mc_frs_errorbars.png'])
        
        figure;
        subplot(2, 1, 1)
        hold on
        plot(inhib_strengths, fr_1(1, :), 'b')
        plot(inhib_strengths, fr_1(2, :), 'r--')
%         plot(inhib_strengths, fr_1(1, :)./(1+inhib_strengths), 'bo-')
%         plot(inhib_strengths, fr_1(2, :)./(1+inhib_strengths), 'ro-')
        hold off
        legend('MC 1', 'MC 2')
        title('Stim 1')
        subplot(2, 1, 2)
        hold on
        plot(inhib_strengths, fr_2(1, :), 'b')
        plot(inhib_strengths, fr_2(2, :), 'r--')
%         plot(inhib_strengths, fr_2(1, :)./(1+inhib_strengths), 'bo-')
%         plot(inhib_strengths, fr_2(2, :)./(1+inhib_strengths), 'ro-')
        hold off
        legend('MC 1', 'MC 2')
        title('Stim 2')
        
        %% process F opt
        % inhib by windows by num_copies by runs
        
        meanF = squeeze(nanmean(squeeze(nanmean(F_opt, 4)), 3));
        mean_pc = squeeze(nanmean(squeeze(nanmean(pc_opt, 4)), 3));
        windowz = repmat(windows, 12, 1);
      
%         meanF_norm = meanF./windowz;
%         
%         figure; 
%         imagesc(windows, inhib_strengths, meanF_norm)
%         title('F opt, normalized by window')
        
        %% F and pc, heat maps
        mf = figure; 
        imagesc(windows, inhib_strengths, meanF)
        title('F opt')
        colorbar;
        xlabel('window size (ms)')
        ylabel('inhib weight')
        hgsave(mf, [target_dir, 'meanF_heat', suffix,'.fig'])
        saveas(mf, [target_dir, 'meanF_heat', suffix,'.png'])
        
        mpc = figure; 
        imagesc(windows, inhib_strengths, mean_pc)
        title('fraction correct')
        xlabel('window size (ms)')
        ylabel('inhib weight')
        colorbar;
        hgsave(mpc, [target_dir, 'meanpc_heat', suffix,'.fig'])
        saveas(mpc, [target_dir, 'meanpc_heat', suffix,'.png'])
        
        %% F and pc vs. inhib. each line is a separate window
        f = figure;
        hold on
        for w_i = 1:length(windows)
            if w_i == 1
                w_1 = plot(inhib_strengths, meanF(:, w_i)./meanF(1, w_i), 'color', window_colorz(w_i, :));
            elseif w_i == length(windows)
                w_end = plot(inhib_strengths, meanF(:, w_i)./meanF(1, w_i), 'color', window_colorz(w_i, :));
            else
                plot(inhib_strengths, meanF(:, w_i)./meanF(1, w_i), 'color', window_colorz(w_i, :));
            end
        end
        hold off
        xlabel('inhibitory weight')
        title('F opt')
        ylabel('proportional change')
%         set(gca, 'Xscale', 'log')
        legend([w_1, w_end], ['bin=', num2str(windows(1))],...
            ['bin=', num2str(windows(end))], 'location', 'northeast')
        hgsave(f, [target_dir, 'meanF', suffix,'.fig'])
        saveas(f, [target_dir, 'meanF', suffix,'.png'])
        
        pcf = figure;
        hold on
        for w_i = 1:length(windows)
            if w_i == 1
                w_1 = plot(inhib_strengths, mean_pc(:, w_i), 'color', window_colorz(w_i, :));
            elseif w_i == length(windows)
                w_end = plot(inhib_strengths, mean_pc(:, w_i), 'color', window_colorz(w_i, :));
            else
                plot(inhib_strengths, mean_pc(:, w_i), 'color', window_colorz(w_i, :));
            end
        end
        hold off
        xlabel('inhibitory weight')
        title('fraction correct')
%         set(gca, 'Xscale', 'log')
        legend([w_1, w_end], ['bin=', num2str(windows(1))],...
            ['bin=', num2str(windows(end))], 'location', 'northeast')
        hgsave(pcf, [target_dir, 'meanpc', suffix,'.fig'])
        saveas(pcf, [target_dir, 'meanpc', suffix,'.png'])
        
        %% F and pc vs. window. each line is a separate inhib.
        fi = figure;
        hold on
        for i_i = 1:length(inhib_strengths)
            if i_i == 1
                i_1 = plot(windows, meanF(i_i, :)./meanF(i_i, 1), 'color', inhib_colorz(i_i, :));
            elseif i_i == length(inhib_strengths)
                i_end = plot(windows, meanF(i_i, :)./meanF(i_i, 1), 'color', inhib_colorz(i_i, :));
            else
                plot(windows, meanF(i_i, :)./meanF(i_i, 1), 'color', inhib_colorz(i_i, :));
            end
        end
        hold off
        xlabel('window size (ms)')
        title('F opt')
        ylabel('proportional change')
%         set(gca, 'Xscale', 'log')
        legend([i_1, i_end], ['w=', num2str(inhib_strengths(1))],...
            ['w=', num2str(inhib_strengths(end))], 'location', 'northwest')
        hgsave(fi, [target_dir, 'meanF_window', suffix,'.fig'])
        saveas(fi, [target_dir, 'meanF_window', suffix,'.png'])
        
        pci = figure;
        hold on
        for i_i = 1:length(inhib_strengths)
            if i_i == 1
                i_1 = plot(windows, mean_pc(i_i, :), 'color', inhib_colorz(i_i, :));
            elseif i_i == length(inhib_strengths)
                i_end = plot(windows, mean_pc(i_i, :), 'color', inhib_colorz(i_i, :));
            else
                plot(windows, mean_pc(i_i, :), 'color', inhib_colorz(i_i, :));
            end
        end
        hold off
        xlabel('window size (ms)')
        title('fraction correct')
%         set(gca, 'Xscale', 'log')
        legend([i_1, i_end], ['w=', num2str(inhib_strengths(1))],...
            ['w=', num2str(inhib_strengths(end))], 'location', 'northwest')
        hgsave(pci, [target_dir, 'meanpc_window', suffix,'.fig'])
        saveas(pci, [target_dir, 'meanpc_window', suffix,'.png'])
       
        
%         f2 = figure;
%         hold on
%         for w_i = 1:length(windows)
%             if w_i == 1
%                 w_1 = plot(inhib_strengths, meanF_norm(:, w_i)./meanF_norm(1, w_i), 'color', window_colorz(w_i, :));
%             elseif w_i == length(windows)
%                 w_end = plot(inhib_strengths, meanF_norm(:, w_i)./meanF_norm(1, w_i), 'color', window_colorz(w_i, :));
%             else
%                 plot(inhib_strengths, meanF_norm(:, w_i)./meanF_norm(1, w_i), 'color', window_colorz(w_i, :));
%             end
%         end
%         hold off
%         xlabel('Inhib. weight')
% %         set(gca, 'Xscale', 'log')
%         legend([w_1, w_end], ['w=', num2str(windows(1))],...
%             ['w=', num2str(windows(end))], 'location', 'northeast')
%         hgsave(f2, [target_dir, 'meanF_norm', suffix,'.fig'])
%         saveas(f2, [target_dir, 'meanF_norm', suffix,'.png'])
        
        %% now we move into the frills section 
        % F opt by mc
        % F opt by mode
        % corrcoefs 
        % covsum_evalues
        % covsum_evecs
        % covsums
        % diffmean 
        
        %% difference in mean by inhibition.
        % similar values for all window sizes.
        f3 = figure;
        hold on
        for w_i = 1:length(windows)
            i1 = plot(inhib_strengths, mean_diffmean(1, :, w_i), 'b');
            plot(inhib_strengths, mean_diffmean(1, :, w_i)./(1+inhib_strengths), 'bo-');
            i2 = plot(inhib_strengths, mean_diffmean(2, :, w_i), 'r');
            plot(inhib_strengths, mean_diffmean(2, :, w_i)./(1+inhib_strengths), 'ro-');
        end
        hold off
        xlabel('inhib')
        ylabel('\Delta \mu')
        legend([i1, i2], 'MC 1', 'MC 2', 'location', 'northeast')
        hgsave(f3, [target_dir, 'diff mu', suffix,'.fig'])
        saveas(f3, [target_dir, 'diff mu', suffix,'.png'])
        
        %% covsums
        % first take the mean, then look at the eigenspectra
        % mean_covsums is Nm x Nm x 12 x 40
        mean_evalues = NaN(Nm, length(inhib_strengths), length(windows));
        mean_evecs = NaN(Nm, Nm, length(inhib_strengths), length(windows));
        % populate these matrices
        for i_i = 1:length(inhib_strengths)
            for w_i = 1:length(windows)
                [ev, ei] = eig(squeeze(mean_covsums(:, :, i_i, w_i)));
                mean_evalues(:, i_i, w_i) = diag(ei);
                mean_evecs(:, :, i_i, w_i) = ev;
            end
        end
        
        %% eigenvalues
        % in all cases for single connections, they are about the same.
        eigf = figure;
        hold on
        for w_i = 1:length(windows)
            if w_i==1
                w_1 = plot(inhib_strengths, mean_evalues(2, :, w_i), 'color', window_colorz(w_i, :));
            elseif w_i==length(windows)
                w_end = plot(inhib_strengths, mean_evalues(2, :, w_i), 'color', window_colorz(w_i, :));
            else 
                plot(inhib_strengths, mean_evalues(2, :, w_i), 'color', window_colorz(w_i, :));
            end
%             plot(inhib_strengths, mean_evalues(1, :, w_i)./(1+inhib_strengths), 'bo-');
%             i2 = plot(inhib_strengths, mean_evalues(2, :, w_i), 'color', window_colorz(w_i, :));
%             plot(inhib_strengths, mean_evalues(2, :, w_i)./(1+inhib_strengths), 'ro-');
        end
        hold off
        xlabel('inhibitory weight')
        ylabel('Eigenvalue')
        legend([w_1, w_end], ['bin=', num2str(windows(1))],...
            ['bin=', num2str(windows(end))], 'location', 'northeast')
        hgsave(eigf, [target_dir, 'eigvalues', suffix,'.fig'])
        saveas(eigf, [target_dir, 'eigvalues', suffix,'.png'])
        
        eigf2 = figure;
        hold on
        for i_i = 1:length(inhib_strengths)
            if i_i==1
                i_1 = plot(windows, squeeze(mean_evalues(1, i_i, :)), 'color', inhib_colorz(i_i, :));
            elseif i_i==length(inhib_strengths)
                i_end = plot(windows, squeeze(mean_evalues(1, i_i, :)), 'color', inhib_colorz(i_i, :));
            else 
                plot(windows, squeeze(mean_evalues(1, i_i, :)), 'color', inhib_colorz(i_i, :));
            end
%             plot(inhib_strengths, mean_evalues(1, :, w_i)./(1+inhib_strengths), 'bo-');
%             i2 = plot(inhib_strengths, mean_evalues(2, :, w_i), 'color', window_colorz(w_i, :));
%             plot(inhib_strengths, mean_evalues(2, :, w_i)./(1+inhib_strengths), 'ro-');
        end
        hold off
        xlabel('window (ms)')
        ylabel('Eigenvalue')
        legend([i_1, i_end], ['w=', num2str(inhib_strengths(1))], ...
            ['w=', num2str(inhib_strengths(end))], 'location', 'northeast')
        hgsave(eigf2, [target_dir, 'eigvalues_window', suffix,'.fig'])
        saveas(eigf2, [target_dir, 'eigvalues_window', suffix,'.png'])
        
        %% time for eigenvectors
        % mean_evecs is Nm x mode x inhib x windows
        
        wi = 40; % chosen window index
        vecf = figure;
        hold on
        for i_i = 1:length(inhib_strengths)
            coord1 = [mean_evecs(1, 1, i_i, wi); -mean_evecs(1, 1, i_i, wi)];
            coord2 = [mean_evecs(2, 1, i_i, wi); -mean_evecs(2, 1, i_i, wi)];
            plot(coord1, coord2, 'color', inhib_colorz(i_i, :), ...
                'linewidth', length(inhib_strengths)-i_i+1)
        end
        xL = xlim;
        yL = ylim;
        line([0 0], yL);  %x-axis
        line(xL, [0 0]);  %y-axis
        hold off
        title('Eigenvector 1')
        hgsave(vecf, [target_dir, 'eigvec1', suffix,'.fig'])
        saveas(vecf, [target_dir, 'eigvec1', suffix,'.png'])
        
        vecf2 = figure;
        hold on
        for i_i = 1:length(inhib_strengths)
            coord1 = [mean_evecs(1, 2, i_i, wi); -mean_evecs(1, 2, i_i, wi)];
            coord2 = [mean_evecs(2, 2, i_i, wi); -mean_evecs(2, 2, i_i, wi)];
            plot(coord1, coord2, 'color', inhib_colorz(i_i, :), ...
                'linewidth', length(inhib_strengths)-i_i+1)
        end
        xL = xlim;
        yL = ylim;
        line([0 0], yL);  %x-axis
        line(xL, [0 0]);  %y-axis
        hold off
        title('Eigenvector 2')
        hgsave(vecf2, [target_dir, 'eigvec2', suffix,'.fig'])
        saveas(vecf2, [target_dir, 'eigvec2', suffix,'.png'])
        
        
    end
end

if ~am_i_on_quest
%% plotting

    poiss_F = squeeze(F_opt(1, :, :));
    netw_F = squeeze(F_opt(2, :, :));

figure;
%     subplot(2, 1, 1)
imagesc(windows, inhib_strengths, netw_F);
title('F opt')
xlabel('window (ms)')
ylabel('inhibition (a.u.)')

figure;
%     subplot(2, 1, 1)
imagesc(windows, inhib_strengths, poiss_F);
title('F opt, Poiss')
xlabel('window (ms)')
ylabel('inhibition (a.u.)')

        poiss_pc = nanmean(squeeze(pc_opt(1, :, :)), 3);
        netw_pc = nanmean(squeeze(pc_opt(2, :, :)), 3);

figure;
%     subplot(2, 1, 1)
imagesc(windows, inhib_strengths, netw_pc);
title('PC')
xlabel('window (ms)')
ylabel('inhibition (a.u.)')

figure;
%     subplot(2, 1, 1)
imagesc(windows, inhib_strengths, poiss_pc);
title('PC, Poiss')
xlabel('window (ms)')
ylabel('inhibition (a.u.)')

%% diff mean
netw_diffmean = squeeze(diff_mean_vector(2, :, :, :));

% across inhibitions
figure;
subplot(2, 1, 1)
imagesc(windows, inhib_strengths, squeeze(netw_diffmean(1, :, :)))
xlabel('window (ms)')
ylabel('inhib')
colorbar;
title('MC 1')
subplot(2, 1, 2)
imagesc(windows, inhib_strengths, squeeze(netw_diffmean(2, :, :)))
xlabel('window (ms)')
ylabel('inhib')
colorbar;
title('MC 2')

covsum_evalues = squeeze(nanmean(covsum_evalues, 6));
covsum_evecs_mean = squeeze(nanmean(covsum_evecs, 7));
netw_evecs = squeeze(covsum_evecs(2, :, :, :, :));

end
