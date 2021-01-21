% 2-MC case as a proof of concept
% to match against the analyical results that when the MCs are laterally
% inhibited by the same GC, that F opt does not appear affected by w

% might be worth calculating mean and covar separately
% to match against the analytical stuff

set(0,'DefaultAxesFontSize',30,'defaultaxeslinewidth',2,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)

T_final = 10000;
inhib_strengths = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 5, 10, 20, 50]; %0:10:50;
osn_scale = 5;

connectivity_type = 'global';
stim_type = 1;

% model parameters
Nm = 2; % number of mitral cells
Ns = 2;
Ng = 1;
N = Ng+Nm;
rp = 0; 
rp_timer = zeros(N, Ns);

dt = 0.02;
times = 0:dt:T_final;
folder_name = '2MC_PoC_norp/';
var_of_interest = inhib_strengths;

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

mc_tau_x = 5; gc_tau_x = 5; osn_tau_x = 1;
tau_x = [mc_tau_x*ones(Nm, Ns); gc_tau_x*ones(Ng, Ns)];

total_runs = length(inhib_strengths);
for tr_i = 1:total_runs
    % do we want to save results to file?
        run_number = dlmread('run_num_network.txt');
        dlmwrite('run_num_network.txt', run_number+1);
        [~, pre_file_location] = save_results(run_number, folder_name);
    file_locations{tr_i} = pre_file_location;
end


weights = [2000, 10000, 10000];
Nt = round(T_final/dt)+1;
background_noise = 1;

if stim_type==1
    stim_widths = 0.3*Nm; 
    stim_center = (1+Nm)/2;
    stim_means = [floor(stim_center-0.1*Nm), ceil(stim_center+0.1*Nm)];
%     stim_means = [floor(stim_center-0.2*Nm), ceil(stim_center+0.2*Nm)];
        
    S = shift_gaussians(Nm, stim_widths, stim_means, background_noise);
end

% n o r m a l i z e 
S= Nm*100*S/sum(S(:)) ;


for vr_i = 1:length(inhib_strengths) %size(variable_run, 2)
    
    inhib_strength = inhib_strengths(vr_i);
    
    %% populate Wmg and Wgm matrices
    wGM = weights(1)*inhib_strength;
    wMG = weights(2); %*inhib_strength;
    
    if strcmp(connectivity_type, 'global')
            Wmg = wMG*ones(Ng, Nm)/Nm;
            Wgm = wGM*ones(Nm, Ng)/Nm;
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
    
        
        file_location = file_locations{vr_i};
        
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
        
        %%%%%%%%%%% START %%%%%%%%%%%
        fprintf('Running......\n');
        
        for count = 1:T_final/dt       
            
            % this gets redefined at each time step
            input = [weights(3)*osn_xbar-... % stimulus to the MCs
                Wgm*xbar(Nm+1:(Nm+Ng),:);... % inhibition
                Wmg*xbar(1:Nm, :)]; % from MCs to GCs
            
            % simulate the OSN spiking
            osn_spike = poissrnd(osn_rate_adjusted + osn_spontaneous_adjusted);

            % these are Nm x Ns
            osn_xbar = osn_xbar+(osn_spike/osn_scale)-dt*osn_xbar/osn_tau_x;
%             osn_xbar_history(:, :, count) = osn_xbar;
            
            % the ones not in refractory get incremented as per usual
            Voltage = Voltage + ...
                dt*(-gL.*(Voltage-EL)+...
                gL*delT*exp((Voltage-VT)/delT)-Hyp+Dep+input)./C;
            
            % the ones where it's nonzero get set to EL
            Voltage(rp_timer~=0)=EL;
            
            spike = Voltage>VT;
            Voltage(spike) = EL; % reset to the resting voltage
%             Hyp(spike) = Hyp(spike) + b;
%             Dep(spike) = Dep(spike) + Isp;
            rp_timer(spike) = rp/dt;
            
            xbar = xbar + spike - (dt*xbar)./tau_x;
            
            Spike_history(:, :, count) = spike; 
%             Voltage_history(:, :, count) = Voltage;
%             
            % this is the element wise maximum
            % prevents negatives 
            rp_timer = max(zeros(N, Ns), rp_timer-1);
        end
        
        fprintf('------END------\n');
        
        Mitral_spike_history = Spike_history(1:Nm, :, 1:end-1);
        Granule_spike_history = Spike_history(Nm+1:(Nm+Ng), :, 1:end-1);
        
        save([file_location, 'data.mat'], 'Mitral_spike_history', ...
            'Granule_spike_history');

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
    copyfile('Fopt_2MC_PoC.m', file_location);
end

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