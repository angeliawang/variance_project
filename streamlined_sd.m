%% useful for Quest:
% bash commands to move folders R15072 through R15142 into new_folder/

% for n in {15072..15142}
% do
%     d="R$n"
%     mv "$d" new_folder/
% done

% dirs = {'single_mix/', 'global_mix/', 'stim_mix/', 'mex_mix/'};
dirs = {'NN_GLOBAL_TEST/'};

am_i_on_quest = 1;
the_frills = 1; % everything beyond Fopt, PC, MC FRs
bool_spectra = 1;
bool_projmean = 1;

set(0, 'DefaultFigureVisible', 'on')
set(0,'DefaultAxesFontSize',30,'defaultaxeslinewidth',2,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)

num_copies = 10; % number of spiking datasets per parameter
inhib_strengths = 0:0.1:1; % for all types of connectivity
windows = 5:5:100;
Nm = 50;
% osn_scales = 5;
var_of_interest = inhib_strengths;
xlabel_of_interest = 'inhibition';

% randomly choose windows.
tessellate_windows = 0;
dt = .02;
T_final = 10000;
Nt = T_final/dt;
itv = 600; % to save this amount 

% number of data points for LDA
% i.e. number of random windows to take
num_samples = 100;

assert(num_samples>= Nm);
% else the covariance matrix is guaranteed singular

% number of times to run this, to smooth things out
num_runs = 100;

for reg = dirs
    
    target_dir = ['Noisy_Neuron/', reg{1}];
    files = dir(target_dir);
    directoryNames = {files([files.isdir]).name};
    directoryNames = directoryNames(~ismember(directoryNames,{'.','..'}));

    % initialize things
    % we use zeros() instead of NaN() because I use a rolling average
    F_opt = zeros(2, length(var_of_interest), length(windows), num_copies*num_runs); % one for Poiss and one for netw
    pc_opt = zeros(size(F_opt));
    mc_frs_1 = zeros(Nm, length(var_of_interest), num_copies);
    mc_frs_2 = zeros(Nm, length(var_of_interest), num_copies);
    
    avg_covsum = zeros(2, Nm, Nm, length(var_of_interest), length(windows));
    
    if bool_spectra
        covsum_evalues = zeros(2, Nm, length(var_of_interest), length(windows), num_copies*num_runs);
    end
    
    if bool_projmean
        diff_means_proj = zeros(2, length(var_of_interest), length(windows), num_copies*num_runs); % one for Poiss and one for netw
    end

    if the_frills
        corrcoefs = zeros(Nm, Nm, length(var_of_interest), length(windows), num_copies*num_runs);
        diff_mean_vector = zeros(2, Nm, length(var_of_interest), length(windows), num_copies*num_runs);
        Fopt_mode = NaN(size(covsum_evalues));
        Fopt_mc = NaN(size(Fopt_mode));
        sample_data = NaN(2, 2, Nm, num_samples, length(var_of_interest), ...
        length(windows));
    end

    for d_i = 1:num_copies*length(var_of_interest) %size(directoryNames, 2)

        d_i
        file_location = [target_dir,directoryNames{d_i}];

        % this will list every .mat file in the directory
    %     files = dir(fullfile(file_location, '*.mat'));

        % use a try catch statement so that if shit is singular, the loop
        % continues onto the next batch of .mat files instead of quitting

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
            
            counter = 1;

                % it's better practice for num_runs to be 1 and the dataset 
                % to be large so as to avoid resampling
                for r_i = 1:num_runs

                data_set = NaN(Nm, Ns, num_samples);

                % for each sample, grab that chunk
                for i = 1:num_samples

                    if tessellate_windows==1
                        % tesselated
                        i1 = (i-1)*(window_steps+1)+1;
                        i2 = (i-1)*(window_steps+1)+1;
                    else
                    % not aligned
                        i1 = randi(num_steps-window_steps);
                        i2 = randi(num_steps-window_steps);
                    end

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
                corrcoefs(:, :, i_i, w_i, (ii_i-1)*num_runs + r_i) = cc;

                % ignore Nms that didn't fire
                channel_keep_mc = find(sum(stim1 + stim2, 2) > 0);
    %             if (length(channel_keep_mc) ~= Nm)
    %                 fprintf(' these channels were kept \n')
    %                 fprintf(' %g ', channel_keep_mc')
    %                 fprintf(' \n')
                    stim1 = stim1(channel_keep_mc, :);
                    stim2 = stim2(channel_keep_mc, :);
    %             end

                if ii_i==1
                    sample_data(2, 1, 1:length(channel_keep_mc), :, i_i, w_i) = ...
                        stim1;
                    sample_data(2, 2, 1:length(channel_keep_mc), :, i_i, w_i) = ...
                        stim2;
                end

                % general poisson spike counts
                % these are Nm by 1
                stim1_mean = mean(stim1, 2);
                stim2_mean = mean(stim2, 2);

                stim1_poiss = poissrnd(repmat(stim1_mean, 1, num_samples));
                stim2_poiss = poissrnd(repmat(stim2_mean, 1, num_samples));

                if ii_i==1
                    sample_data(1, 1, 1:length(channel_keep_mc), :, i_i, w_i) = ...
                        stim1_poiss;
                    sample_data(1, 2, 1:length(channel_keep_mc), :, i_i, w_i) = ...
                        stim2_poiss;
                end
    %             stim1_poiss = NaN(size(stim1));
    %             stim2_poiss = NaN(size(stim2)); 
    %             for ns_i = 1:num_samples
    %                 stim1_poiss(:, ns_i) = poissrnd(stim1_mean);
    %                 stim2_poiss(:, ns_i) = poissrnd(stim2_mean);
    %             end

                % poiss
                diff_means_poiss = mean(stim1_poiss, 2)-mean(stim2_poiss, 2);
                cov_poiss = diag(mean(stim1_poiss, 2)+mean(stim2_poiss, 2));
                w_poiss = cov_poiss\diff_means_poiss;
                c_poiss = dot(w_poiss, (stim1_mean+stim2_mean))/2;
                P1c_poiss = w_poiss'*stim1_poiss;
                P2c_poiss = w_poiss'*stim2_poiss;
                pc_opt(1, i_i, w_i, (ii_i-1)*num_runs + r_i) = ...
                    sum([P1c_poiss>c_poiss P2c_poiss<c_poiss])/num_samples/2;
                            % for poisson, F opt scales linearly with wdow size, so
                % there is a normalizing factor
                F_opt(1, i_i, w_i, (ii_i-1)*num_runs + r_i) = (dot(w_poiss, diff_means_poiss))^2/(w_poiss'*cov_poiss*w_poiss)/wdow;

                if bool_spectra
                    [vecs, vals] = eig(cov_poiss);
                    [vals_s, ind] = sort(diag(vals));
                    vecs = vecs(:, ind);
                    covsum_evalues(1, 1:length(channel_keep_mc), i_i, w_i, (ii_i-1)*num_runs + r_i) = vals_s/wdow;
%                     current = covsum_evecs(1, :, :, i_i, w_i);
%                     covsum_evecs(1, 1:length(channel_keep_mc), 1:length(channel_keep_mc), i_i, w_i) = ((counter-1)*current+vecs)/counter;
                end

                if bool_projmean
                    diff_means_proj(1, i_i, w_i, (ii_i-1)*num_runs + r_i) = dot(diff_means_poiss,w_poiss);
                end
                
                if the_frills
                    avg_covsum(1, :, :, i_i, w_i) = (squeeze(avg_covsum(1, :, :, i_i, w_i))*(counter-1) + cov_poiss)/counter;

                    for mode_i = 1:length(channel_keep_mc)
                        mc_i = channel_keep_mc(mode_i);
                        Fopt_mode(1, mode_i, i_i, w_i, (ii_i-1)*num_runs + r_i) = ...
                            dot(vecs(:, mode_i), diff_means_poiss)^2/vals_s(mode_i);
                        Fopt_mc(1, mc_i, i_i, w_i, (ii_i-1)*num_runs + r_i) = ...
                            diff_means_poiss(mode_i)*w_poiss(mode_i);
                    end

                    diff_mean_vector(1, 1:length(channel_keep_mc), i_i, w_i, (ii_i-1)*num_runs + r_i) = diff_means_poiss;
                end

                % netw
                diff_means = mean(stim1, 2)-mean(stim2, 2);
                cov_sum = cov(stim1')+cov(stim2');
                w = cov_sum\diff_means;
                c = dot(w, (stim1_mean+stim2_mean))/2;
                P1c = w'*stim1;
                P2c = w'*stim2;
                pc_opt(2, i_i, w_i, (ii_i-1)*num_runs + r_i) = ...
                    sum([P1c>c P2c<c])/num_samples/2;
                F_opt(2, i_i, w_i, (ii_i-1)*num_runs + r_i) = (dot(w, diff_means))^2/(w'*cov_sum*w)/wdow;

                if bool_spectra
                    [vecs, vals] = eig(cov_sum);
                    [vals_s, ind] = sort(diag(vals));
                    vecs = vecs(:, ind);
                    covsum_evalues(2, 1:length(channel_keep_mc), i_i, w_i, (ii_i-1)*num_runs + r_i) = vals_s/wdow;
                end
                
                if bool_projmean
                    diff_means_proj(2, i_i, w_i, (ii_i-1)*num_runs + r_i) = dot(diff_means, w);
                end

                if the_frills
                    avg_covsum(2, :, :, i_i, w_i) = (squeeze(avg_covsum(2, :, :, i_i, w_i))*(counter-1) + cov_sum)/counter;
                    counter = counter + 1;
                    
                    for mode_i = 1:length(channel_keep_mc)
                        mc_i = channel_keep_mc(mode_i);
                        Fopt_mode(2, mode_i, i_i, w_i, (ii_i-1)*num_runs + r_i) = ...
                            dot(vecs(:, mode_i), diff_means)^2/vals_s(mode_i);
                        Fopt_mc(2, mc_i, i_i, w_i, (ii_i-1)*num_runs + r_i) = ...
                            diff_means(mode_i)*w(mode_i);
                    end

                    diff_mean_vector(2, 1:length(channel_keep_mc), i_i, w_i, (ii_i-1)*num_runs + r_i) = diff_means;

                end

                end

    %             ni = length(var_of_interest);
    %             nw = length(windows);
    %             if ii_i==1
    %                 % stim1 is Nm by num_samples
    %                 figure(101)
    %                 hold all
    %                 subplot(nw, ni, (w_i-1)*ni+i_i)
    %                 scatter(stim1(14, :), stim1(34, :), 20, 'r', 'filled')
    %                 hold on
    %                 scatter(stim2(14, :), stim2(34, :), 10, 'b', 'filled')
    %                 if w_i==3
    %                     ylim([0 25])
    %                     xlim([0 25])
    %                 else
    %                     ylim([0 10])
    %                     xlim([0 10])
    %                 end
    % %                 axis square
    % %                 title(['netw, ii ', num2str(i_i), ' wdow ', num2str(wdow)])
    % 
    %                 % poiss scatterpoints not filled.
    %                 figure(102)
    %                 subplot(nw, ni, (w_i-1)*ni+i_i)
    %                 scatter(stim1_poiss(14, :), stim1_poiss(34, :), 20, 'r')
    %                 hold on
    %                 scatter(stim2_poiss(14, :), stim2_poiss(34, :), 10, 'b')
    %                 if w_i==3
    %                     ylim([0 25])
    %                     xlim([0 25])
    %                 else
    %                     ylim([0 10])
    %                     xlim([0 10])
    %                 end
    % %                 axis square
    % %                 title(['poiss, ii ', num2str(i_i), ' wdow ', num2str(wdow)])
    %             end
        end
        catch
        continue
        end

    end

	%covsum_evecs = squeeze(nanmean(covsum_evecs, 6));
    % corrcoefs is Nm x Nm x inhibs x windows x num_copies
    % covsum_evalues is Nm x inhibs x windows x num_copies

    save(strcat(target_dir, 'LDA.mat'), 'F_opt', 'pc_opt', 'mc_frs_1', ...
        'mc_frs_2')
    if bool_spectra
        save(strcat(target_dir, 'Spectra.mat'), 'covsum_evalues')
    end
    if bool_projmean
        save(strcat(target_dir, 'projmean.mat'), 'diff_means_proj')
    end
    if the_frills
        save(strcat(target_dir, 'sample_data.mat'), 'sample_data')
        save(strcat(target_dir, 'Frills.mat'), 'corrcoefs', 'avg_covsum', 'diff_mean_vector', 'diff_means_proj', ...
            'covsum_evalues', 'Fopt_mode', 'Fopt_mc')
    end

end

