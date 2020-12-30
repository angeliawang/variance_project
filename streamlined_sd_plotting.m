% regs = {'single_mix/', 'global_mix/', 'stim_mix/', 'mex_mix/', ...
%     'single_skew/', 'global_skew/', 'stim_skew/', 'mex_skew/'};
regs = {'single_mix/'};

set(0, 'DefaultFigureVisible', 'on')
set(0,'DefaultAxesFontSize',30,'defaultaxeslinewidth',2,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)

suffix = '';
universal_colorbar = 0;
bool_Fopt = 1;
bool_pc = 1;
bool_mcfr = 1;
bool_windowdep = 1;
bool_cov_evalues = 1;
bool_meanvec = 1;
bool_Fopt_mc = 1;
% if applicable
    ii = 4; % chosen inhibitory index
    first_window = 1;
    wi = 10; % chosen window index
        
num_copies = 50; % number of spiking datasets per parameter
inhib_strengths = 0:100:1000; % for all types of connectivity
% windows = 1:2:20;
windows = 5:5:200;
Nm = 50;
% osn_scales = 5;
var_of_interest = inhib_strengths;
xlabel_of_interest = 'inhibition';
num_inhibs = length(inhib_strengths);
num_windows = length(windows);

w_colorz = varycolor(num_windows);
i_colorz = varycolor(num_inhibs);
m_colorz = varycolor(Nm);

for reg = regs
    
    target_dir = ['AdExIF/', reg{1}];
    
    if bool_Fopt
        load([target_dir, 'LDA', suffix, '.mat']);
        F_opt(abs(F_opt)>30) = NaN;
    end
    
    if bool_cov_evalues
        frills = load([target_dir, 'Frills', suffix, '.mat']);
        covsum_evalues = frills.covsum_evalues;
        covsum_evalues(covsum_evalues==0) = NaN;
        mean_evalues = squeeze(mean(covsum_evalues, 5));
        netw_evalues = squeeze(mean_evalues(2, :, ii, :));
        poiss_evalues = squeeze(mean_evalues(1, :, ii, :));
    end
    if bool_meanvec
%         load([target_dir, 'diffmeans', suffix, '.mat']);
        diff_mean_vector = frills.diff_mean_vector;
        mean_meanvec = squeeze(mean(diff_mean_vector, 5));
        netw_meanvec = squeeze(mean_meanvec(2, :, :, :));
        poiss_meanvec = squeeze(mean_meanvec(1, :, :, :));
    end
    if bool_Fopt_mc
        cov_evecs = frills.covsum_evecs;
        
        % 50x50x11x50: Nm by mode by inhibs by sample
        selected_evecs = squeeze(cov_evecs(2, :, :, :, wi, :));
        selected_evals = squeeze(covsum_evalues(2, :, :, wi, :));
        selected_diffmean = squeeze(diff_mean_vector(2, :, :, wi, :));
    end
    
%     F_chunk = F_opt(:, :, 4:end, :);
%     max(F_chunk(:))
%     min(F_chunk(:))
%     continue

    % load([target_dir, 'sample_data.mat']);
    % sample_mean = squeeze(mean(sample_data, 4));
    % sample_std = squeeze(std(sample_data, 0, 4));

    % these are indices of the corresponding vector.
%     selected_windows = [1, 3, 40];
    selected_inhibs = [3, 11];

     if bool_Fopt
    %% F opt and PC
    poiss_F = nanmean(squeeze(F_opt(1, :, :, :)), 3);
    netw_F = nanmean(squeeze(F_opt(2, :, :, :)), 3);
    std_poiss_F = nanstd(squeeze(F_opt(1, :, :, :)), 0, 3);
    std_netw_F = nanstd(squeeze(F_opt(2, :, :, :)), 0, 3);

        f = figure;
    %     subplot(2, 1, 1)
        imagesc(windows(first_window:end), inhib_strengths, netw_F(:, first_window:end));
    %     title([target_dir, 'F opt netw'], 'Interpreter', 'none')
        title('F opt')
        xlabel('window (ms)')
        ylabel('inhibition (a.u.)')
        colorbar;
    %     caxis([0, 0.5])
        hgsave(f, [target_dir, 'F_netw_heatmap', suffix, '.fig'])
        saveas(f, [target_dir, 'F_netw_heatmap', suffix, '.png'])
    %     if universal_colorbar
    %         caxis([0, 0.29])
    %     end
        f2 = figure;
        imagesc(windows(first_window:end), inhib_strengths, poiss_F(:, first_window:end));
        title('F opt, Poiss')
        xlabel('window (ms)')
        ylabel('inhibition (a.u.)')
        colorbar;
    %     caxis([0, 0.5])
        hgsave(f2, [target_dir, 'F_poiss_heatmap', suffix, '.fig'])
        saveas(f2, [target_dir, 'F_poiss_heatmap', suffix, '.png'])
    end

    if bool_pc
        poiss_pc = nanmean(squeeze(pc_opt(1, :, :, :)), 3);
        netw_pc = nanmean(squeeze(pc_opt(2, :, :, :)), 3);

        f3 = figure;
    %     subplot(2, 1, 1)
        imagesc(windows(first_window:end), inhib_strengths, netw_pc(:, first_window:end));
        title('PC')
        xlabel('window (ms)')
        ylabel('inhibition (a.u.)')
        colorbar;
        caxis([0.5, 1])
        hgsave(f3, [target_dir, 'pc_netw_heatmap', suffix, '.fig'])
        saveas(f3, [target_dir, 'pc_netw_heatmap', suffix, '.png'])

        f4 = figure;
        imagesc(windows(first_window:end), inhib_strengths, poiss_pc(:, first_window:end));
        title('PC, Poiss')
        xlabel('window (ms)')
        ylabel('inhibition (a.u.)')
        colorbar;
        caxis([0.5, 1])
        hgsave(f4, [target_dir, 'pc_poiss_heatmap', suffix, '.fig'])
        saveas(f4, [target_dir, 'pc_poiss_heatmap', suffix, '.png'])
    end

        %% MC FRs
        if bool_mcfr
%  
            mc_fr_m1 = mean(mc_frs_1, 3);
            mc_fr_m1_std = std(mc_frs_1, 0, 3);
            mc_fr_m2 = mean(mc_frs_2, 3);
            mc_fr_m2_std = std(mc_frs_2, 0, 3);
% 
            ColorSet = varycolor(length(var_of_interest));
            f5 = figure;
            hold on
            for i_i = 1:length(var_of_interest)
            plot(mc_fr_m1(:, i_i), 'color', ColorSet(i_i, :))
            plot(mc_fr_m2(:, i_i), '--', 'color', ColorSet(i_i, :))
            end
            ylabel('kHz')
            xlabel('MC index')
            title('MC FR')
            xlim([0 51])
            hgsave(f5, [target_dir, 'mc_frs_colored.fig'])
            saveas(f5, [target_dir, 'mc_frs_colored.png'])
            % matlab2tikz('figurehandle',f2,'filename',[target_dir, 'mc_frs_colored.tex'] ,'standalone', true);

            % another visual of MC FRs
%             ColorSet = varycolor(length(var_of_interest));
%             f2 = figure;
%             hold on
%             for i_i = 1:length(var_of_interest)
%             errorbar(mc_fr_m1(:, i_i), mc_fr_m1_std(:, i_i), 'color', ColorSet(i_i, :))
%             errorbar(mc_fr_m2(:, i_i), mc_fr_m2_std(:, i_i), '--', 'color', ColorSet(i_i, :))
%             end
%             ylabel('kHz')
%             xlabel('MC index')
%             title([target_dir, 'MC FR'], 'Interpreter', 'none')
%             hgsave([target_dir, 'mc_frs_colored_error.fig'])
%             saveas(f2, [target_dir, 'mc_frs_colored_error.epsc'])
            % matlab2tikz('figurehandle',f2,'filename',[target_dir, 'mc_frs_colored_error.tex'] ,'standalone', true);
        end
        
    if bool_windowdep
        F_chunk_poiss = squeeze(poiss_F(ii, first_window:end));
        F_chunk_netw = squeeze(netw_F(ii, first_window:end));
        F_std_poiss = squeeze(std_poiss_F(ii, first_window:end));
        F_std_netw = squeeze(std_netw_F(ii, first_window:end));

        f6 = figure;
        hold on
        errorbar(windows(first_window:end), F_chunk_poiss, F_std_poiss, 'bo-')
        errorbar(windows(first_window:end), F_chunk_netw, F_std_netw, 'ro-')
        hold off
        legend('Poiss', 'Netw', 'location', 'northwest')
        xlabel('window (ms)')
        ylabel('F opt')
        hgsave(f6, [target_dir, 'Fopt_by_window.fig'])
        saveas(f6, [target_dir, 'Fopt_by_window.png'])
    end
    
    if bool_cov_evalues
        
        f7 = figure;
        hold on
        for w_i = first_window:num_windows
            plot(netw_evalues(:, w_i), 'o-', 'Color', w_colorz(w_i-first_window+1, :))
        end
        xlim([0, 50])
        hold off
        title('\Sigma eigenvalues')
        xlabel('mode')
        ylabel('magnitude')
        hgsave(f7, [target_dir, 'cov_evals_netw', suffix, '.fig'])
        saveas(f7, [target_dir, 'cov_evals_netw', suffix, '.png'])
        
        f8 = figure;
        hold on
        for w_i = first_window:num_windows
            plot(poiss_evalues(:, w_i), 'o-', 'Color', w_colorz(w_i-first_window+1, :))
        end
        xlim([0, 50])
        hold off
        title('Poiss \Sigma eigenvalues')
        xlabel('mode')
        ylabel('magnitude')
        hgsave(f8, [target_dir, 'cov_evals_poiss', suffix, '.fig'])
        saveas(f8, [target_dir, 'cov_evals_poiss', suffix, '.png'])
        
        f13 = figure;
        hold on
        for i_i = 1:num_inhibs
            plot(squeeze(mean_evalues(2, :, i_i, wi)), 'Color', i_colorz(i_i, :));
        end
        xlim([0, 50])
        hold off
        title('\Sigma eigenvalues, by inhib')
        xlabel('mode')
        ylabel('magnitude')
        hgsave(f13, [target_dir, 'cov_evals_inhib', suffix, '.fig'])
        saveas(f13, [target_dir, 'cov_evals_inhib', suffix, '.png'])
        
        f14 = figure;
        hold on
        for i_i = 1:num_inhibs
            plot(squeeze(mean_evalues(1, :, i_i, wi)), 'Color', i_colorz(i_i, :));
        end
        xlim([0, 50])
        hold off
        title('\Sigma eigenvalues Poiss, by inhib')
        xlabel('mode')
        ylabel('magnitude')
        hgsave(f8, [target_dir, 'cov_evals_inhib_poiss', suffix, '.fig'])
        saveas(f8, [target_dir, 'cov_evals_inhib_poiss', suffix, '.png'])
        
        f15 = figure;
        hold on
        for i_i = 1:num_inhibs
            plot(1./squeeze(mean_evalues(2, :, i_i, wi)), 'Color', i_colorz(i_i, :));
        end
        xlim([0, 50])
        hold off
        title('$\Sigma \frac{1}{\lambda}$, by inhib', 'interpreter', 'latex')
        xlabel('mode')
        ylabel('magnitude')
        hgsave(f15, [target_dir, 'cov_evals_inhib_recip', suffix, '.fig'])
        saveas(f15, [target_dir, 'cov_evals_inhib_recip', suffix, '.png'])
        
        f16 = figure;
        hold on
        for i_i = 1:num_inhibs
            plot(1./squeeze(mean_evalues(1, :, i_i, wi)), 'Color', i_colorz(i_i, :));
        end
        xlim([0, 50])
        hold off
        title('$\Sigma \frac{1}{\lambda}$ Poiss, by inhib', 'interpreter', 'latex')
        xlabel('mode')
        ylabel('magnitude')
        hgsave(f16, [target_dir, 'cov_evals_inhib_recip_poiss', suffix, '.fig'])
        saveas(f16, [target_dir, 'cov_evals_inhib_recip_poiss', suffix, '.png'])
        
    end
    
    if bool_meanvec        
        
        f9 = figure;
        hold on
        for w_i = first_window:num_windows
            plot(netw_meanvec(:, ii, w_i), 'o-', 'Color', w_colorz(w_i-first_window+1, :))
        end
        xlim([0, 50])
        hold off
        title('\Delta \mu, by window')
        xlabel('MC')
        ylabel('magnitude')
        hgsave(f9, [target_dir, 'netw_meanvec_bywindow', suffix, '.fig'])
        saveas(f9, [target_dir, 'netw_meanvec_bywindow', suffix, '.png'])
        
        f10 = figure;
        hold on
        for w_i = first_window:num_windows
            plot(poiss_meanvec(:, ii, w_i), 'o-', 'Color', w_colorz(w_i-first_window+1, :))
        end
        xlim([0, 50])
        hold off
        title('\Delta \mu, Poiss, by window')
        xlabel('MC')
        ylabel('magnitude')
        hgsave(f10, [target_dir, 'poiss_meanvec_bywindow', suffix, '.fig'])
        saveas(f10, [target_dir, 'poiss_meanvec_bywindow', suffix, '.png'])
        
        f11 = figure;
        hold on
        for i_i = 1:11
            plot(netw_meanvec(:, i_i, wi), 'o-', 'Color', i_colorz(i_i, :))
        end
        xlim([0, 50])
        hold off
        title('\Delta \mu, by inhib')
        xlabel('MC')
        ylabel('magnitude')
        hgsave(f11, [target_dir, 'netw_meanvec_byinhib', suffix, '.fig'])
        saveas(f11, [target_dir, 'netw_meanvec_byinhib', suffix, '.png'])
        
        f12 = figure;
        hold on
        for i_i = 1:11
            plot(poiss_meanvec(:, i_i, wi), 'o-', 'Color', i_colorz(i_i, :))
        end
        xlim([0, 50])
        hold off
        title('\Delta \mu Poiss, by inhib')
        xlabel('MC')
        ylabel('magnitude')
        hgsave(f12, [target_dir, 'poiss_meanvec_byinhib', suffix, '.fig'])
        saveas(f12, [target_dir, 'poiss_meanvec_byinhib', suffix, '.png'])
    end
    
    if bool_Fopt_mc
        % Nm by inhib by samples by modes
        Fopt_components = zeros(50, 11, 50, 50);
        
        % this is very inefficient
        for mc_i = 1:Nm
            for inhib_i = 1:num_inhibs
                for sample_i = 1:50
                    for mode_i = 1:Nm
                    Fopt_components(mc_i, inhib_i, sample_i, mode_i) = ...
                        (selected_diffmean(mc_i, inhib_i, sample_i)*...
                        selected_evecs(mc_i, mode_i, inhib_i, sample_i))^2/...
                        selected_evals(mode_i, inhib_i, sample_i);
                    end
                end
            end
        end
        
        mean_Fopt_components = squeeze(mean(Fopt_components, 4));
        Fopt_by_mc = squeeze(sum(mean_Fopt_components, 3)); % sum over modes

        % line plots
        figure;
        hold on
        for m_i = 1:Nm
            plot(Fopt_by_mc(m_i, :)/windows(wi)^2, 'color', m_colorz(m_i, :))
        end
        hold off
        
        % heat map, which tbh I prefer
        figure;
        imagesc(inhib_strengths, 1:Nm, Fopt_by_mc/windows(wi)^2)
        colorbar;
        ylabel('MC')
        xlabel('inhibition (a.u.)')
        title('Contribution of each MC to F opt')
        
        % sanity check, sum over MCs and we should get the whole F opt
        % the stored F opt divides by the window. also, eigenvalues are
        % div'd by the window
        figure;
        plot(inhib_strengths, sum(Fopt_by_mc, 1)/windows(wi)^2)
        xlabel('inhib strengths')
    end
        
    %% diff means as a vector
    % % 2 x Nm x inhibs x windows x num_copies*num_runs
    % diff_mean_vector = squeeze(mean(diff_mean_vector, 5));
    % diff_mean_netw = squeeze(diff_mean_vector(2, :, :, :));
    % diff_mean_poiss = squeeze(diff_mean_vector(1, :, :, :));
    % 
    % colors = varycolor(length(inhib_strengths));
    % 
    % for w_i = selected_windows %1:length(windows)
    %     wdow = windows(w_i);
    % 
    %     figure
    %     subplot(2, 1, 1)
    %     hold on
    %     for i = 1:length(inhib_strengths)
    %        plot(diff_mean_netw(:, i, w_i)./wdow, 'color', colors(i, :))
    %     end
    %     hold off
    %     xlabel('MC')
    %     title(['diff means, wdow', num2str(wdow)])
    %     subplot(2, 1, 2)
    %     hold on
    %     for i = 1:length(inhib_strengths)
    %        plot(diff_mean_poiss(:, i, w_i)./wdow, 'color', colors(i, :))
    %     end
    %     hold off
    %     xlabel('MC')
    %     hgsave([target_dir, num2str(wdow), '_diffmeans.fig'])
    % end
    %     
    %     
    %     %% diff means as a projection
    %     diff_means_proj = nanmean(diff_means_proj, 4);
    %     % this is just F opt itself, but not scaled
    % 
    %     for w_i = selected_windows %1:length(windows)
    %         wdow = windows(w_i);
    % 
    %         figure;
    %         hold on
    %         plot(inhib_strengths, diff_means_proj(2, :, w_i), 'b', 'linewidth', 3)
    %         plot(inhib_strengths, diff_means_proj(1, :, w_i), 'r')
    %         hold off
    %         xlabel('inhibition')
    %         ylabel('projection of diff means')
    %         legend('netw', 'poiss')
    %         title(['diff means projection, wdow', num2str(wdow)])
    %         hgsave([target_dir, num2str(wdow), '_diffmeansproj.fig'])
    %     end
    % 
    % %% eigvalues of the cov mtx
    % 
    % covsum_evalues = squeeze(nanmean(covsum_evalues, 5));
    % 
    % cov_evalues_netw = squeeze(covsum_evalues(2, :, :, :));
    % cov_evalues_poiss = squeeze(covsum_evalues(1, :, :, :));
    % 
    % for w_i = selected_windows %1:length(windows)
    %     wdow = windows(w_i);
    % 
    %     figure;
    %     subplot(2, 1, 1)
    %     imagesc(inhib_strengths, 1:Nm, log(cov_evalues_netw(:, :, w_i)))
    %     colorbar;
    %     ylabel('mode')
    %     xlabel('inhibition')
    %     title(['cov eigvalues, wdow', num2str(wdow)])
    %     subplot(2, 1, 2)
    %     imagesc(inhib_strengths, 1:Nm, log(cov_evalues_poiss(:, :, w_i)))
    %     colorbar;
    %     ylabel('mode')
    %     xlabel('inhibition')
    %     hold off
    %     title(['poiss cov eigvalues, wdow', num2str(wdow)])
    %     hgsave([target_dir, num2str(wdow), '_cov_eigvalues.fig'])
    % 
    % end
    % 
    %     %% contributions to F opt by mode and by MC
    %     Fopt_mode(isinf(Fopt_mode))=NaN;
    %     Fopt_mode(abs(Fopt_mode)>100) = NaN;
    %     Fopt_mode = squeeze(nanmean(Fopt_mode, 5));
    %     Fopt_mode_netw = squeeze(Fopt_mode(2, :, :, :));
    %     Fopt_mode_poiss = squeeze(Fopt_mode(1, :, :, :));
    %     
    %     Fopt_mc(isinf(Fopt_mc))=NaN;
    %     Fopt_mc(abs(Fopt_mc)>100) = NaN;
    %     Fopt_mc = squeeze(nanmean(Fopt_mc, 5));
    %     Fopt_mc_netw = squeeze(Fopt_mc(2, :, :, :));
    %     Fopt_mc_poiss = squeeze(Fopt_mc(1, :, :, :));
    % 
    %     for w_i = selected_windows %1:length(windows)
    %         wdow = windows(w_i);
    % 
    %         figure;
    %         subplot(2, 1, 1)
    %         imagesc(inhib_strengths, 1:Nm, squeeze(Fopt_mode_netw(:, :, w_i)))
    %         colorbar;
    %         xlabel('inhibition')
    %         ylabel('mode')
    %         title(['mode contributions to Fopt, wdow', num2str(wdow)])
    %         subplot(2, 1, 2)
    %         imagesc(inhib_strengths, 1:Nm, squeeze(Fopt_mode_poiss(:, :, w_i)))
    %         colorbar;
    %         xlabel('inhibition')
    %         ylabel('mode')
    %         hgsave([target_dir, num2str(wdow), '_Fopt_mode_map.fig'])
    %         
    %     end
    %         
    %     for w_i = selected_windows %1:length(windows)
    %         wdow = windows(w_i);
    %         figure;
    %         subplot(2, 1, 1)
    %         imagesc(inhib_strengths, 1:Nm, squeeze(Fopt_mc_netw(:, :, w_i)))
    %         colorbar;
    %         ylabel('mc')
    %         xlabel('inhibition')
    %         title(['mc contributions to Fopt, wdow', num2str(wdow)])
    %         subplot(2, 1, 2)
    %         imagesc(inhib_strengths, 1:Nm, squeeze(Fopt_mc_poiss(:, :, w_i)))
    %         colorbar;
    %         ylabel('mc')
    %         xlabel('inhibition')
    %         hgsave([target_dir, num2str(wdow), '_Fopt_mc_map.fig'])
    %         
    %     end
    % 
    % %% eigenvectors of the cov mtx
    % % 2 x Nm x Nm x num_inhibs x windows
    % 
    % covsum_evecs = squeeze(nanmean(covsum_evecs, 6));
    % covsum_evecs_poiss = squeeze(covsum_evecs(1, :, :, :, :));
    % covsum_evecs_netw = squeeze(covsum_evecs(2, :, :, :, :));
    % % these are Nm x Nm x inhibs x windows
    %     
    % % let's start with ordered by mode (ascending corresponding eigenvalue)
    % colors = varycolor(Nm);
    % for w_i = selected_windows
    %     wdow = windows(w_i);
    % %     figure
    % %     for i_i = 1:length(inhib_strengths)
    % %     subplot(4, 3, i_i)
    % %     dis = squeeze(covsum_evecs_netw(:, :, i_i, w_i));
    % %     imagesc(dis)
    % %     ylabel('mode')
    % %     title(['inhib', num2str(inhib_strengths(i_i))])
    % %     colorbar;
    % %     end
    %     % not great
    %     
    % %     figure
    % %     for i_i = 1:length(inhib_strengths)
    % %     subplot(4, 3, i_i)
    % %     dis = squeeze(covsum_evecs_netw(:, :, i_i, w_i));
    % %     hold on
    % %     for bop = 1:size(dis, 1)
    % %         plot(dis(:, bop), 'color', colors(bop, :))
    % %     end
    % %     hold off
    % %     sgtitle(['Sorted by ascending eigenvalue, wdow', num2str(wdow)])
    % %     title(['inhib', num2str(inhib_strengths(i_i))])
    % %     end
    %     
    %     % sort by contribution 
    %     figure
    %     for i_i = 1:length(inhib_strengths)
    %         [~, s_ind] = sort(Fopt_mode(2, :, i_i, w_i));
    %         
    %         subplot(4, 3, i_i)
    %         dis = squeeze(covsum_evecs_netw(:, s_ind, i_i, w_i));
    %         hold on
    %         for bop = 1:3 %size(dis, 1)
    %             plot(dis(:, bop), 'color', colors(bop, :))
    %         end
    %         hold off
    %         title(['inhib', num2str(inhib_strengths(i_i))])
    %     end
    %     sgtitle(['Sorted by contribution, wdow', num2str(wdow)])
    %     hgsave([target_dir, num2str(wdow), '_evecs_sorted.fig'])
    % end
    % 
    %% sample mean at selected inhibs
    % for i_i = selected_inhibs
    %     inhib = inhib_strengths(i_i);
    %     
    %     figure;
    %     % netw stim 1
    %     subplot(2, 1, 1)
    %     imagesc(squeeze(sample_mean(2, 1, :, i_i, :)))
    %     colorbar;
    %     ylabel('Nm')
    %     xlabel('window')
    %     title('spike avg, 100 samples')
    %     subplot(2, 1, 2)
    %     imagesc(squeeze(sample_mean(1, 1, :, i_i, :)))
    %     colorbar;
    %     ylabel('Nm')
    %     xlabel('window')
    %     title('Poiss')
    %     hgsave([target_dir, num2str(inhib), '_samplespikes.fig'])
    % end

end

%% a bit of resampling business of questionable use
% w_i = 3;
% data1 = squeeze(sample_data(1, 1, :, :, :, w_i));
% data2 = squeeze(sample_data(1, 2, :, :, :, w_i));
% wdows = 5:5:20;
% wdow = wdows(w_i);
% % these should be Nm x 2000samples x inhibs
% num_inhibs = size(data1, 3);
% 
% num_samples = 1000:1000:8000;
% num_trials = 50;
% 
% Fopt_by_sample = NaN(length(num_samples), num_inhibs, num_trials);
% 
% for nt_i = 1:num_trials
% 
% for ns_i = 1:length(num_samples)
%     num_samps = num_samples(ns_i);
% 
%     samps1 = data1(:, randi(8000, 1, num_samps), :);
%     samps2 = data2(:, randi(8000, 1, num_samps), :);
%     % these are Nm x num_samples x inhibs
%     
%     for i_i = 1:num_inhibs
%         sample1 = squeeze(samps1(:, :, i_i));
%         sample2 = squeeze(samps2(:, :, i_i));
%         %Nm x num_samples
%         
%         diff_means = mean(sample1, 2)-mean(sample2, 2);
%         cov_sum = cov(sample1') + cov(sample2');     
%         w = cov_sum\diff_means;
%         
%         % divide by window size 5
%         Fopt_by_sample(ns_i, i_i, nt_i) = (dot(w, diff_means))^2/(w'*cov_sum*w)/wdow;
%         
%     end
% end
% end
% 
% mean_Fopt = squeeze(mean(Fopt_by_sample, 3));
% std_Fopt = squeeze(std(Fopt_by_sample, 0, 3));
% 
% colors = varycolor(length(inhib_strengths));
% 
% figure;
% hold on
% for i_i = 3
% errorbar(num_samples, mean_Fopt(:, i_i), std_Fopt(:, i_i), 'color', colors(i_i, :))
% end
% hold off
% 
% for i_i = 3
%     figure;
%     hold on
%     errorbar(num_samples, mean_Fopt(:, i_i)-0.033, std_Fopt(:, i_i), 'color', colors(i_i, :))
%     plot(num_samples, .005224*(num_samples/1000).^(-.5))
%     hold off
%     title(['inhibition', num2str(inhib_strengths(i_i))])
%     set(gca, 'YScale', 'log')
%     set(gca, 'XScale', 'log')
% end
%     
% figure;
% imagesc(num_samples, inhib_strengths, mean_Fopt')
% colorbar;
% ylabel('inhibition')
% xlabel('number of samples')
% title('single skew window 10 Poiss, F opt by sample')

