dirs = {'2MC_PoC/'};

set(0, 'DefaultFigureVisible', 'on')
set(0,'DefaultAxesFontSize',30,'defaultaxeslinewidth',2,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)

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
num_samples = 50;

runs = 5;

assert(num_samples>= Nm);
% else the covariance matrix is guaranteed singular

% number of times to run this, to smooth things out

for reg = dirs
        
    target_dir = ['AdExIF/', reg{1}];
    files = dir(target_dir);
    directoryNames = {files([files.isdir]).name};
    directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

    % usually these are the same
    % there is the odd run where quest timed out and we didn't hit all the
    % inhib weights
    sup = min(length(var_of_interest), length(directoryNames)/num_copies);
        
    % initialize things
    % we use zeros() instead of NaN() because I use a rolling average
    F_opt = zeros(2, sup, length(windows), num_copies, runs); % one for Poiss and one for netw
    pc_opt = zeros(size(F_opt));
    mc_frs_1 = zeros(Nm, sup);
    mc_frs_2 = zeros(Nm, sup);
    
    covsums = NaN(Nm, Nm, sup, length(windows), num_copies, runs);
    covsum_evalues = zeros(2, Nm, sup, length(windows), num_copies, runs);
    covsum_evecs = NaN(2, Nm, Nm, sup, length(windows), num_copies, runs);

        corrcoefs = zeros(Nm, Nm, sup, length(windows), num_copies, runs);
        diff_mean_vector = zeros(2, Nm, sup, length(windows), num_copies, runs);
        Fopt_mode = NaN(size(covsum_evalues));
        Fopt_mc = NaN(size(Fopt_mode));
        
    for d_i = 1:size(directoryNames, 2) 

        d_i
        file_location = [target_dir,directoryNames{d_i}];

        full_path = [file_location, '/data.mat'];
    %     full_path = [file_location, '/', file.name];
        data_mc = load(full_path);
        
        Mitral_spike_history = data_mc.Mitral_spike_history;
        Nm = size(Mitral_spike_history, 1);
        Ns = size(Mitral_spike_history, 2);
        num_steps = size(Mitral_spike_history, 3);

        % Do some stuff

        % inhib_strength indices
        i_i = ceil(d_i/num_copies);
        ii_i = mod(d_i, num_copies); % from 1 to num_copies
        if ii_i==0
            ii_i = num_copies;
        end

        % mc_frs are Nm by num_inhibs by num_copies
        % and contain *FIRING RATES*
        mc_frs_1(:, i_i, ii_i) = sum(Mitral_spike_history(:, 1, :), 3)/(T_final); 
        mc_frs_2(:, i_i, ii_i) = sum(Mitral_spike_history(:, 2, :), 3)/(T_final); 

        try 
        %% for each wdow
        for w_i = 1:length(windows)
            wdow = windows(w_i);
            window_steps = round(wdow/dt)-1;

                data_set = NaN(Nm, Ns, num_samples);
                
                for r_i = 1:runs

                % for each sample, grab that chunk
                for i = 1:num_samples

                    i1 = randi(num_steps-window_steps);
                    i2 = randi(num_steps-window_steps);

                    % these contain SPIKE COUNTS in a designated wdow
                    data_set(:, 1, i) = sum(Mitral_spike_history(:, 1, i1:i1+window_steps), 3);
                    data_set(:, 2, i) = sum(Mitral_spike_history(:, 2, i2:i2+window_steps), 3);
                end

                % these are Nm by num_samples
                % and contain spike totals
                stim1 = squeeze(data_set(:, 1, :));
                stim2 = squeeze(data_set(:, 2, :));

                cc = corrcoef(stim1') + corrcoef(stim2');
                % corrcoefs is Nm x Nm x inhibs x windows x
                % (num_copies*num_runs)
                corrcoefs(:, :, i_i, w_i, (ii_i-1)*runs + r_i) = cc;

                % ignore Nms that didn't fire
                channel_keep_mc = find(sum(stim1 + stim2, 2) > 0);
    %             if (length(channel_keep_mc) ~= Nm)
    %                 fprintf(' these channels were kept \n')
    %                 fprintf(' %g ', channel_keep_mc')
    %                 fprintf(' \n')
                    stim1 = stim1(channel_keep_mc, :);
                    stim2 = stim2(channel_keep_mc, :);
    %             end

%                 if ii_i==1
%                     sample_data(2, 1, 1:length(channel_keep_mc), :, i_i, w_i) = ...
%                         stim1;
%                     sample_data(2, 2, 1:length(channel_keep_mc), :, i_i, w_i) = ...
%                         stim2;
%                 end
                
                % these are Nm by 1
                stim1_mean = mean(stim1, 2);
                stim2_mean = mean(stim2, 2);
                
                % generate poisson spike counts
                stim1_poiss = poissrnd(repmat(stim1_mean, 1, num_samples));
                stim2_poiss = poissrnd(repmat(stim2_mean, 1, num_samples));

%                 if ii_i==1 && r_i==1
%                     sample_data(1, 1, 1:length(channel_keep_mc), :, i_i, w_i) = ...
%                         stim1_poiss;
%                     sample_data(1, 2, 1:length(channel_keep_mc), :, i_i, w_i) = ...
%                         stim2_poiss;
%                 end
    %             stim1_poiss = NaN(size(stim1));
    %             stim2_poiss = NaN(size(stim2)); 
    %             for ns_i = 1:num_samples
    %                 stim1_poiss(:, ns_i) = poissrnd(stim1_mean);
    %                 stim2_poiss(:, ns_i) = poissrnd(stim2_mean);
    %             end

                % now that samples are generated, let's convert them to FRs
                % in kHz
                stim1 = stim1/wdow;
                stim2 = stim2/wdow;
                stim1_poiss = stim1_poiss/wdow;
                stim2_poiss = stim2_poiss/wdow;
                
                % poiss
                diff_means_poiss = mean(stim1_poiss, 2)-mean(stim2_poiss, 2);
                cov_poiss = diag(mean(stim1_poiss, 2)+mean(stim2_poiss, 2));
                w_poiss = cov_poiss\diff_means_poiss;
                c_poiss = dot(w_poiss, (stim1_mean+stim2_mean))/2;
                P1c_poiss = w_poiss'*stim1_poiss;
                P2c_poiss = w_poiss'*stim2_poiss;
                pc_opt(1, i_i, w_i, (ii_i-1)*runs + r_i) = ...
                    sum([P1c_poiss>c_poiss P2c_poiss<c_poiss])/num_samples/2;
                            % for poisson, F opt scales linearly with wdow size, so
                % there is a normalizing factor
%                 F_opt(1, i_i, w_i) = (dot(w_poiss, diff_means_poiss))^2/(w_poiss'*cov_poiss*w_poiss)/wdow;
                F_opt(1, i_i, w_i, (ii_i-1)*runs + r_i) = (dot(w_poiss, diff_means_poiss))^2/(w_poiss'*cov_poiss*w_poiss);

                    [vecs, vals] = eig(cov_poiss);
                    [vals_s, ind] = sort(diag(vals));
                    vecs = vecs(:, ind);
                    covsum_evalues(1, 1:length(channel_keep_mc), i_i, w_i, (ii_i-1)*runs + r_i) = vals_s/wdow;
                    covsum_evecs(1, 1:length(channel_keep_mc), ...
                        1:length(channel_keep_mc), i_i, w_i, (ii_i-1)*runs + r_i) = vecs;
                
                    for mode_i = 1:length(channel_keep_mc)
                        mc_i = channel_keep_mc(mode_i);
                        Fopt_mode(1, mode_i, i_i, w_i, (ii_i-1)*runs + r_i) = ...
                            dot(vecs(:, mode_i), diff_means_poiss)^2/vals_s(mode_i);
                        Fopt_mc(1, mc_i, i_i, w_i, (ii_i-1)*runs + r_i) = ...
                            diff_means_poiss(mode_i)*w_poiss(mode_i);
                    end

                    diff_mean_vector(1, 1:length(channel_keep_mc), ...
                        i_i, w_i, (ii_i-1)*runs + r_i) = diff_means_poiss;

                % netw
                diff_means = mean(stim1, 2)-mean(stim2, 2);
                cov_sum = cov(stim1')+cov(stim2');
                w = cov_sum\diff_means;
                c = dot(w, (stim1_mean+stim2_mean))/2;
                P1c = w'*stim1;
                P2c = w'*stim2;
                pc_opt(2, i_i, w_i, (ii_i-1)*runs + r_i) = ...
                    sum([P1c>c P2c<c])/num_samples/2;
%                 F_opt(2, i_i, w_i) = (dot(w, diff_means))^2/(w'*cov_sum*w)/wdow;
                F_opt(2, i_i, w_i, (ii_i-1)*runs + r_i) = (dot(w, diff_means))^2/(w'*cov_sum*w);
                
                covsums(1:length(channel_keep_mc), 1:length(channel_keep_mc), ...
                    i_i, w_i, (ii_i-1)*runs + r_i) = cov_sum;

                [vecs, vals] = eig(cov_sum);
                [vals_s, ind] = sort(diag(vals));
                vecs = vecs(:, ind);
                covsum_evalues(2, 1:length(channel_keep_mc), i_i, w_i, ...
                    (ii_i-1)*runs + r_i) = vals_s/wdow;
                covsum_evecs(2, 1:length(channel_keep_mc), ...
                    1:length(channel_keep_mc), i_i, w_i, (ii_i-1)*runs + r_i) = vecs;
                
                for mode_i = 1:length(channel_keep_mc)
                    mc_i = channel_keep_mc(mode_i);
                    Fopt_mode(2, mode_i, i_i, w_i, (ii_i-1)*runs + r_i) = ...
                        dot(vecs(:, mode_i), diff_means)^2/vals_s(mode_i);
                    Fopt_mc(2, mc_i, i_i, w_i, (ii_i-1)*runs + r_i) = ...
                        diff_means(mode_i)*w(mode_i);
                end

                diff_mean_vector(2, 1:length(channel_keep_mc), ...
                        i_i, w_i, (ii_i-1)*runs + r_i) = diff_means;
            end
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

end

if ~am_i_on_quest
%% plotting
colorz = varycolor(length(inhib_strengths));

% this is a one time thing to address an indexing issue I had
% F_opt_meaned = NaN(2, 7, 40);
% for type_i = [1, 2]
%     for w_i = 1:40
%         for i_i = 1:7
%             start_i = (i_i-1)*50 + 1;
%             end_i = i_i*50;
%             F_opt_meaned(type_i, i_i, w_i) = ...
%                 nanmean(F_opt(type_i, start_i:end_i, w_i));
%         end
%     end
% end

figure;
hold on
for i = 1:length(inhib_strengths)
plot(mc_frs_1(:, i), 'color', colorz(i, :))
plot(mc_frs_2(:, i), '--', 'color', colorz(i, :))
end
xlim([0.5, 2.5])
hold off

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