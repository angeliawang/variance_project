function spiking_network_rp()

% whereby the "_rp" indicates refractory period 

% this bullcrap is apparently necessary before calling parfor
% delete(gcp('nocreate'))
% parpool(2)

set(0,'DefaultAxesFontSize',30,'defaultaxeslinewidth',2,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)

%% THINGS WE CHANGE FOR TROUBLESHOOTING
runs = 50; % number of copies per parameter
T_final = 10000; % in milliseconds
inhib_strengths = [0, 1, 5, 10, 20, 50, 75]; %0:50:1000; %[0, 0.4]; %0:50:1000; %:0.1:1;
osn_scales = 5; %1:10; % times the freq and reduce height of xbar
scale_both = 0; % if excit strength should scale too

connectivity_type = 'global'; % global % stim % mex
stim_type = 4; % 1 = shifted gauss, 2 = skewed gauss, % 4 = mixtures
if stim_type==4
    mix_stretch = 0.2; %0:0.1:0.5; % dist from center, default is 0.15
end
sniff = 0;
weight_noise = 1; % want noise in weight matrices?

am_i_on_quest = 1;

Nm = 50; % number of mitral cells
Ns = 2; % number of stimuli
GM_ratio = 1; % GCs per MC
Ng = GM_ratio*Nm;
N = Ng + Nm; % total # of cells

% in milliseconds
rp = 5; 
rp_timer = zeros(N, Ns);

dt = 0.02;
times = 0:dt:T_final;
folder_name = 'global_mix_smaller_inhibs/';
var_of_interest = inhib_strengths;
osn_scale = osn_scales(1);

% 0 to 1. 
% 0 means weight matrix is just the mean. 1 is fully stretched.
weight_stretch = 0.5; 

% set(0, 'DefaultFigureVisible', 'on')
% note: this will not suppress the figures made by separate functions
%       specifically, shift_gaussians, plot_now, and plot_history

%% model parameters
% from Table 1 in Clopath et al. 2010
C = 281; % membrane capacitance (pF)
gL = 30; EL = -70.6; % leak conductance (nS) and resting potential (mV).
delT = 2; % "slope factor" (mV)
Vthresh = -50.4; % threshold at rest (mV)
tau_ad = 144; % adaptation time constant (ms)
a = 4; b = 0.805; % subthresh adaptation and spike triggered adaptation
% Isp = 400; % spike current after a spike (pA)
Isp = 40;
tau_z = 40; % spike current time constant (ms)
tau_T = 50; % threshold potential time constant (ms)
Vmax = 30.4; % threshold potential after a spike (mV)
    % I've found that qualitative results from snowbird are sort of hard to
    % replicate with an adaptive firing threshold.
% Vmin = -110;

% these are from visual cortex experiments
% section b, Table 1 in Clopath et al. 2010 
mc_tau_x = 5; gc_tau_x = 5; osn_tau_x = 1;
tau_x = [mc_tau_x*ones(Nm, Ns); gc_tau_x*ones(Nm, Ns)];

if sniff
%% experimental parameters
sniff_freq = 2; % Hz
sniff_period = 1000/sniff_freq; % milliseconds
end

% bool_fft = 0;

%% variables of interest that affect the spiking data
% variables_run = {inhib_strengths, ... % inhibitory strengths
%                  [-90], ...  % Vmin, default is -90mV
%                  [stim_type]};    % types of stims:
%                           1: shifted Gaussians
%                           2: skew Gaussians with the same primacy set

% columns are lists of combinations of variables
% rows correspond to types of variable
% variable_run = repackage(variables_run);

%% how do you want your sniff cycle in the morning?

% towards the goal of having a more realistically shaped sniff cycle
% model as the difference of two exponentials
if sniff
    times_mod = mod(times, sniff_period);
    tau_1 = 1*sniff_period/10; tau_2 = 2*sniff_period/10; % in milliseconds
    sniff_intensity = exp(-times_mod/tau_2) - exp(-times_mod/tau_1);
    sniff_intensity = sniff_intensity/max(sniff_intensity); % normalize to 1
else
    % if you instead want a step function sniff cycle (off and on)
    % sniff_intensity = 1-heaviside(times_mod-50);

    % if you just want a constant stimulus, i.e. no sniff cycle
    sniff_intensity = ones(size(times));
end

total_runs = runs*length(var_of_interest);

% this shite is to take care of overlap in directories
% so instead, we make all the directories at the start before anything gets
% run
file_locations = cell(total_runs, 1);
saveresults = 1;
for tr_i = 1:total_runs
    % do we want to save results to file?
    if saveresults
        run_number = dlmread('run_num_network.txt');
        dlmwrite('run_num_network.txt', run_number+1);
        [~, pre_file_location] = save_results(run_number, folder_name);
    end
    file_locations{tr_i} = pre_file_location;
end

% GC-|MC, MC->GC, OSN->MC
weights = [10000, 5000, 10000000];
Nt = round(T_final/dt)+1;

background_noise = 0;

if stim_type==1
    stim_widths = 0.3*Nm; 
    stim_center = (1+Nm)/2;
    stim_means = [floor(stim_center-0.1*Nm), ceil(stim_center+0.1*Nm)];
%     stim_means = [floor(stim_center-0.2*Nm), ceil(stim_center+0.2*Nm)];
        
    S = shift_gaussians(Nm, stim_widths, stim_means, background_noise);
elseif stim_type==2 % skew Gaussians
    S = skew_gaussians(Nm, 1);
elseif stim_type==3 % stim block diags
    brick_size = Nm/5;
    curvy_boi = skew_gaussians(Nm, background_noise);
    S = blockify(curvy_boi, brick_size);
elseif stim_type==4 % mixtures
    stim_widths = 0.2*Nm; 
    stim_center = (1+Nm)/2;
    stim_means = [floor(stim_center-0.2*Nm), ceil(stim_center+0.2*Nm)];
    
    S0 = shift_gaussians(Nm, stim_widths, stim_means, background_noise);
    
    ratio1 = 1-mix_stretch;
    ratio2 = 1+mix_stretch;
    S = [(ratio1*S0(:, 1) + ratio2*S0(:, 2))/2, (ratio2*S0(:, 1) + ratio1*S0(:, 2))/2];

end

% n o r m a l i z e 
S = Nm*100*S/sum(S(:));

for vr_i = length(var_of_interest) %size(variable_run, 2)
    
%     osn_scale = var_of_interest(vr_i);
    inhib_strength = var_of_interest(vr_i);
    
    %% populate Wmg and Wgm matrices
    wGM = weights(1)*inhib_strength;
    wMG = weights(2); %*inhib_strength;
    if scale_both
        wMG = wMG/inhib_strength;
    end
    
    % as a sanity check, it's good to make sure the weights are what you
    % think you wrote them to be.
%     figure;
%     imagesc(Wmg);
    
    % S, desired firing rates
    Ns = size(S, 2);
    osn_rate = S*osn_scale; % in Hz
    osn_spontaneous = 0; % in Hz %%%%%%%
    osn_rate_adjusted = dt*osn_rate/1000; % expected rate per time step
    osn_spontaneous_adjusted = osn_spontaneous*dt/1000;
    
    for r_i = 1:runs % per parameter set
        % can submit a factor of 24 nodes(or processors)
        
        fprintf(strcat('Run', num2str(r_i), '\n'))
        
	file_location = file_locations{(vr_i-1)*runs+r_i};
	if mod(r_i, 2)==1
		stream = RandStream('mt19937ar', 'Seed', (r_i+1)/2);
        	RandStream.setGlobalStream(stream);
	
        %% NOISE
%         spont_rate = 10; %Hz
%         stim_latencies = 20*(100 - floor(S+5*randn(size(S)))); % countssss
%         block_sniff = zeros(Nm, Ns, 40/dt);
%         for nmi = 1:Nm
%             block_sniff(nmi, 1, stim_latencies(nmi, 1)) = 1;
%             block_sniff(nmi, 2, stim_latencies(nmi, 2)) = 1;
%         end

        % if there is noise in the weight matrix, this needs to go inside
        % the run loop.
        if strcmp(connectivity_type, 'stim_block')
            % then have the MCs --> GCs block-diagonally
            brick = ones(Nm/5); % heheheh
            Wmg = wMG*blkdiag(brick, brick, brick, brick, brick);
            Wgm = wGM*eye(Nm)*Nm;
        elseif strcmp(connectivity_type, 'stim_blockdiag')
            brick = ones(Nm/5); % heheheh
            base_Wmg = wMG*(weight_S*weight_S'/3);
            Wmg = base_Wmg.*blkdiag(brick, brick, brick, brick, brick);

            % pentadiag for Wgm
            Wgm = wGM*(eye(Nm)+diag(ones(Nm-1, 1), 1)+...
                diag(ones(Nm-1, 1), -1)+...
                diag(ones(Nm-2, 1), 2)+...
                diag(ones(Nm-2, 1), -2))*Nm/2;

            % now we stretch
            Wmg_diff = Wmg-mean(Wmg(:));
            Wmg = mean(Wmg(:))+weight_stretch*Wmg_diff;
        elseif strcmp(connectivity_type, 'stim')
            % for stim-based weights
            % scale appropriately and self-mulitply to form the weights
            scale_S = S-min(S(:));
            weight_S = scale_S(:, 1)+scale_S(:, 2);  
            weight_S = weight_S/(max(weight_S(:)));

            Wmg = wMG*(weight_S*weight_S');

            % pentadiag
            Wgm = wGM*(eye(Nm)+diag(ones(Nm-1, 1), 1)+...
                diag(ones(Nm-1, 1), -1)+...
                diag(ones(Nm-2, 1), 2)+...
                diag(ones(Nm-2, 1), -2));

            % now we stretch
            Wmg_diff = Wmg-mean(Wmg(:));
            Wmg = mean(Wmg(:))+weight_stretch*Wmg_diff;

            Wmg = Wmg/Nm;
            Wgm = Wgm/5;
        elseif strcmp(connectivity_type,'custom')
            wmm = load('wmm.mat');
            wmm = wmm.Wmm;
            Wmg = (wmm/max(wmm(:)) + eye(Nm))*wMG;

            % pentadiag
            Wgm = wGM*(eye(Nm)+diag(ones(Nm-1, 1), 1)+...
                diag(ones(Nm-1, 1), -1)+...
                diag(ones(Nm-2, 1), 2)+...
                diag(ones(Nm-2, 1), -2))*Nm/2;
        elseif strcmp(connectivity_type, 'global')
            Wmg = wMG*ones(Ng, Nm)/Nm;
            Wgm = wGM*ones(Nm, Ng)/Nm;
        elseif strcmp(connectivity_type, 'single')
            % assumes Nm=Ng;
            Wmg = wMG*eye(Nm);
            Wgm = wGM*eye(Ng);
        elseif strcmp(connectivity_type, 'mex')
            % excit
            % pentadiag
            Wmg = wMG*(eye(Nm)+diag(ones(Nm-1, 1), 1)+...
                diag(ones(Nm-1, 1), -1)+...
                diag(ones(Nm-2, 1), 2)+...
                diag(ones(Nm-2, 1), -2));
            
            % inhib
            amp_inh = wGM/2;
            amp_exc = wGM;
            width_exc = Nm/25;
            width_inh = Nm/5;
            Wgm = zeros(Nm, Ng);
            for im = 1:Nm
                for ig = 1:Ng
                    Wgm(im, ig) = -amp_exc*exp(-(ig-im)^2/width_exc^2)+...
                        amp_inh*exp(-(ig-im)^2/width_inh^2);
                end
            end
            % don't forget that Wgm gets SUBTRACTED
        end

        if weight_noise
            Wmg = Wmg.*(rand(Nm)>0.1);
            Wgm = Wgm.*(rand(Nm)>0.1);
        end
       
	end
 
        %% initialize things
        % THINGS YOU NEED REGARDLESS
        % but initialize them as zeros instead of NaNs
        % so that we can skip ahead when there are RPs
        Spike_history = zeros(N, Ns, Nt); % spikes for both stims at all times
        mc_input_history = zeros(Nm, Ns, Nt);
        gc_input_history = zeros(Ng, Ns, Nt);
%         osn_xbar_history = NaN(Nm, Ns, Nt);
%         Voltage_history = NaN(N, Ns, Nt);
        
        % resets everything
        Voltage = EL*ones(N, Ns);
        xbar = zeros(N, Ns);
        osn_xbar = zeros(Nm, Ns);
        Hyp = zeros(size(Voltage)); % w_ad, hyperpolarization current
        Dep = zeros(size(Voltage)); % after-spike depolarization current
        VT = Vthresh*ones(size(Voltage));
 
	rp_timer = zeros(N, Ns);
       
        %%%%%%%%%%% START %%%%%%%%%%%
        fprintf('Running......\n');
        
        for count = 1:T_final/dt       
            
            % this gets redefined at each time step
            input = [weights(3)*osn_xbar-... % stimulus to the MCs
                Wgm*xbar(Nm+1:(Nm+Ng),:);... % inhibition
                Wmg*xbar(1:Nm, :)]; % from MCs to GCs
            
%             input = [S*weights(3)-... % stimulus to the MCs
%                 Wgm*xbar(Nm+1:(Nm+Ng),:);... % inhibition
%                 Wmg*xbar(1:Nm, :)]; % from MCs to GCs
            
            if (~am_i_on_quest)
                mc_input_history(:, :, count) = input(1:Nm, :);
                gc_input_history(:, :, count) = input(Nm+1:N, :);
            end
            
            % simulate the OSN spiking
            osn_spike = poissrnd(osn_rate_adjusted*sniff_intensity(count) + osn_spontaneous_adjusted);
            
%             stim_index = mod(count, 40/dt);
%             if stim_index ==0
%                 stim_index = 40/dt;
%             end
%             noise = rand(Nm, Ns)<(spont_rate*dt/1000);
%             osn_spike = block_sniff(:, :, stim_index) + noise;
            
            % these are Nm x Ns
            osn_xbar = osn_xbar+(osn_spike/osn_scale)-dt*osn_xbar/osn_tau_x;
%             osn_xbar_history(:, :, count) = osn_xbar;
            
            % the ones not in refractory get incremented as per usual
            Voltage = Voltage + ...
                dt*(-gL.*(Voltage-EL)+...
                gL*delT*exp((Voltage-VT)/delT)-Hyp+Dep+input)./C;
            
            % the ones where it's nonzero get set to EL
            Voltage(rp_timer~=0)=EL;
%             VT = VT + dt*(Vthresh-VT)/tau_T;
%             Voltage(Voltage<Vmin) = Vmin;
            
            spike = Voltage>VT;
            Voltage(spike) = EL; % reset to the resting voltage
%             Hyp(spike) = Hyp(spike) + b;
%             Dep(spike) = Dep(spike) + Isp;
            rp_timer(spike) = rp/dt;
            
            xbar = xbar + spike - (dt*xbar)./tau_x;
            
            Spike_history(:, :, count) = spike; 
%             Voltage_history(:, :, count) = Voltage;
%             
%             if mod(count, time_til_peak/dt)==0
%                 % resets everything
%                 Voltage = EL*ones(N, 2);
%                 xbar = zeros(N, 2);
%                 osn_xbar = zeros(Nm, Ns);
% %                 VT = Vthresh*ones(size(Voltage));
%             end
            
            % this is the element wise maximum
            % prevents negatives 
            rp_timer = max(zeros(N, Ns), rp_timer-1);
        end
        
        fprintf('------END------\n');
        
        Mitral_spike_history = Spike_history(1:Nm, :, 1:end-1);
        Granule_spike_history = Spike_history(Nm+1:(Nm+Ng), :, 1:end-1);
%         osn_xbar_history = osn_xbar_history(:, :, 1:end-1);
%         Voltage_history = Voltage_history(:, :, 1:end-1);
        
        if (~am_i_on_quest)
            
        spikess = sum(Mitral_spike_history, 3);
        figure;
        hold on;
        plot(spikess(:, 1))
        plot(spikess(:, 2))
        hold off
        
        % if you want raster plots to be plotted.
        latency_hists;
%         mc_frs = sum(Mitral_spike_history, 3)/(T_final);
%         figure;
%         hold on
%         plot(mc_frs(:, 1))
%         plot(mc_frs(:, 2))
%         ylabel('kHz')
%         hold off

        % sum of spikes div by Nm, Ns, and total time in milliseconds
        avg_mc_fr = sum(Mitral_spike_history(:))/(Ns*Nm*T_final)
        avg_gc_fr = sum(Granule_spike_history(:))/(Ns*Ng*T_final)
        
%         figure;
%         imagesc(squeeze(osn_xbar_history(:, 1, :)))
%         title('osn_xbar, stim 1')
        
        mc_input_1 = squeeze(mc_input_history(:, 1, :));
        figure;
        imagesc(mc_input_1)
        avg_mc_input_1 = mean(mc_input_1(:));
        title(['MC stim1, avg input=', num2str(avg_mc_input_1)])
%         
%         figure;
%         imagesc(squeeze(mc_input_history(:, 2, :)))
%         title('input to MCs, stim 2')
%         
        gc_input_1 = squeeze(gc_input_history(:, 1, :));
        figure;
        imagesc(gc_input_1)
        avg_gc_input_1 = mean(gc_input_1(:));
        title(['GC stim1, avg input=', num2str(avg_gc_input_1)])
%         
%         figure;
%         imagesc(squeeze(Voltage_history(1:Nm, 1, :)));
%         title('mc voltage history')
%         
%         figure;
%         imagesc(squeeze(Voltage_history(Nm+1:end, 1, :)));
%         title('gc voltage history')

        end
       
%         
        parsave([file_location, 'data.mat'], Mitral_spike_history, ...
            Granule_spike_history);
%         save([file_location, 'osn_xbar_history.mat'], 'osn_xbar_history');
%         figure; imagesc(squeeze(osn_xbar_history(:, 1, :)));
    
        %% Spectra
        % we only do this for network-driven spiking, and only once
        % does not depend on current window size
%         SS1 = squeeze(Mitral_spike_history(:, 1, :));  % for spike_spectra.m
%         if bool_fft && r_i == 1
%             windows = variable_diagnostic(2, :);
%             spike_spectra(SS1, dt, inhib_strength, file_location, windows) %separate script
%         end
    end

end

if ~am_i_on_quest
    mc_fr = squeeze(sum(Mitral_spike_history, 3))/T_final;

    % another visual of MC FRs
    figure;
    hold on
    plot(mc_fr(:, 1), 'b')
    plot(mc_fr(:, 2), 'b--')
    ylabel('kHz')
    xlabel('MC index')
    title('MC FR')
end

end

function parsave(save_here, Mitral_spike_history, Granule_spike_history)
    save(save_here, 'Mitral_spike_history', 'Granule_spike_history');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = shift_gaussians(Nm, stim_stddevs, stim_means, background_noise)
% assumes the following exist:   
%           stim_stdevs, a scalar
%           stim_means, a length-2 vector containing desired means
%           background noise, exactly what it sounds like
% outputs:  S, an NMx2 matrix where each column is an odor
%           and an imagesc of S

x = 1:Nm;

% make some Gaussi bois
S = NaN(Nm, length(stim_means));
S(:, 1) = exp(-(x-stim_means(1)).^2/stim_stddevs^2); % odor 1
S(:, 2) = exp(-(x-stim_means(2)).^2/stim_stddevs^2); % odor 2

% add some noise
S = S + background_noise;

end

function S = skew_gaussians(Nm, background_noise)

x = 1:Nm;
center = (Nm+1)/2;

gaussian = @(x) (1/sqrt((2*pi))*exp(-x.^2/2));
skewedgaussian = @(x,alpha) 2*gaussian(x).*normcdf(alpha*x);

% shift the center rightward
S = NaN(Nm, 2);
stretch = 30;
temp1 =  skewedgaussian((x-center)/stretch, 5);
temp2 = skewedgaussian(-(x-center)/stretch, 5);
[~, ind1] = max(temp1);
[~, ind2] = max(temp2);
S(:, 1) = skewedgaussian((x+(ind1-center)-center)/stretch, 5);
S(:, 2) = skewedgaussian(-(x-(center-ind2)-center)/stretch, 5);

% add some noise and normalized if desired
S = S + background_noise;

end

function S = blockify(curvy_S, brick_size)

S = NaN(size(curvy_S));
Nm = size(curvy_S, 1);
checkpoints = 1:brick_size:Nm;
for c_i = 1:length(checkpoints)
    start_i = checkpoints(c_i);
    end_i = checkpoints(c_i)+brick_size-1;
    S(start_i:end_i, 1) = min(curvy_S(start_i:end_i, 1));
    S(start_i:end_i, 2) = min(curvy_S(start_i:end_i, 2));
end

end

function S = steps(Nm)
S = NaN(Nm, 2);

S(:, 1) = [0.2*ones(0.2*Nm, 1); 2*ones(0.6*Nm, 1); zeros(0.2*Nm, 1)];
S(:, 2) = [zeros(0.2*Nm, 1); 2*ones(0.6*Nm, 1);0.2*ones(0.2*Nm, 1)];
end

function [file_tag, file_location] = save_results(run_number, folder)

% different code runs if date_hour_in is given
%     date_hour= clock;
%     date_hour= strcat(num2str(date_hour(1),'%02i'),'_',...
%         num2str(date_hour(2),'%02i'),'_',...
%         num2str(date_hour(3),'%02i'),'_',...
%         num2str(date_hour(4),'%02i'),'_',...
%         num2str(date_hour(5),'%02i'));
    file_tag = strcat('R',num2str(run_number)); %,'_',date_hour);
 
fprintf('This is run '); fprintf('%s',file_tag); fprintf('\n');

% set up subdirectories
mkdir(['AdExIF/', folder, file_tag]);
file_location = ['AdExIF/', folder, file_tag,'/'];
    copyfile('spiking_network_rp.m', file_location);
end

function coordinates = repackage(vars)
% inputs:       vars, a cell array of arbitrary length of 1D double arrays
%               contents not necessarily of the same length
% outputs:      coordinates, a matrix whose columns are coordinates of 
%               every combination of the indices represented in vars
%  e.g. repackage({[1, 2, 3], [4, 5]}) should return 
%      1     2     3     1     2     3
%      4     4     4     5     5     5
% e.g. repackage({[1], [2, 3], [4, 5]}) should return
%      1     1     1     1
%      2     3     2     3
%      4     4     5     5
% e.g. repackage({[1]}) should return {[1]}

%% AFTER STRESSING ABOUT THIS FOR A WHILE, I DISCOVERED THAT MATLAB
% HAS A BUILTIN FUNCTION THAT DOES EXACTLY WHAT I WANT

num_vars = length(vars);
coordinates = vars{1};

if num_vars > 1
    for v_i = 2:num_vars
        coordinates = combvec(coordinates, vars{v_i});
    end
end

end
