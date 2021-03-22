% 2-MC LINEAR RATE MODEL

set(0,'DefaultAxesFontSize',30,'defaultaxeslinewidth',2,...
    'defaultlinelinewidth',2.,'defaultpatchlinewidth',1.5)

num_copies = 20;
T_final = 20000;
inhib_strengths = [0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 5, 10, 20, 50];
%[0, 0.1, 0.2, 0.4, 0.6, 0.8, 1, 2, 5, 10, 20, 50]; %0:10:50;

connectivity_type = 'global';
stim_type = 1;

% model parameters
Nm = 2; % number of mitral cells
Ns = 2;
Ng = 2;
N = Ng+Nm;

dt = 0.02;
times = 0:dt:T_final;
folder_name = '2MC_global_linear_exploremore/';
var_of_interest = inhib_strengths;

mc_tau_x = 1; gc_tau_x = 1; %dt/10;
tau_x = [mc_tau_x*ones(Nm, Ns); gc_tau_x*ones(Ng, Ns)];

total_runs = length(inhib_strengths)*num_copies;
file_locations = cell(total_runs, 1);
for tr_i = 1:total_runs
    % do we want to save results to file?
        run_number = dlmread('run_num_network.txt');
        dlmwrite('run_num_network.txt', run_number+1);
        [~, pre_file_location] = save_results(run_number, folder_name);
    file_locations{tr_i} = pre_file_location;
end

Nt = round(T_final/dt)+1;
background_noise = 1;

if stim_type==1
    stim_widths = 0.3*Nm; 
    stim_center = (1+Nm)/2;
    stim_means = [floor(stim_center-0.1*Nm), ceil(stim_center+0.1*Nm)];
     stim_means = [floor(stim_center-0.2*Nm), ceil(stim_center+0.2*Nm)];
        
    S = shift_gaussians(Nm, stim_widths, stim_means, background_noise);

% n o r m a l i z e 
S= S/sum(S(:));
elseif stim_type==2
	S = [65.3130, 34.6870; 65.3130, 34.6870];
end

weights = [1, 1];
for vr_i = 1:length(inhib_strengths) %size(variable_run, 2)
    
    inhib_strength = inhib_strengths(vr_i);
    
    %% populate Wmg and Wgm matrices
    wGM = weights(1)*inhib_strength;
    wMG = weights(2); %*inhib_strength;
    
    if strcmp(connectivity_type, 'global')
            Wmg = wMG*ones(Ng, Nm)/Nm;
            Wgm = wGM*ones(Nm, Ng)/Nm;
    elseif strcmp(connectivity_type, 'single')
        Wmg = wMG*eye(Ng);
        Wgm = wGM*eye(Ng);
    end
    
    % as a sanity check, it's good to make sure the weights are what you
    % think you wrote them to be.
%     figure;
%     imagesc(Wmg);
    
    % S, desired firing rates
   Ns = 2;
 
       for r_i = 1:num_copies 
        file_location = file_locations{(vr_i-1)*num_copies + r_i};
        
        %% initialize things
        % THINGS YOU NEED REGARDLESS
        % but initialize them as zeros instead of NaNs
        % so that we can skip ahead when there are RPs
        Voltage_history = zeros(N, Ns, Nt); % spikes for both stims at all times
%         osn_xbar_history = NaN(Nm, Ns, Nt);
%         Voltage_history = NaN(N, Ns, Nt);
        
        % resets everything
        Voltage = zeros(N, Ns);
        
        %%%%%%%%%%% START %%%%%%%%%%%
        fprintf('Running......\n');
        
        for count = 1:T_final/dt       
            
            % this gets redefined at each time step
            input = [S-Wgm*Voltage(Nm+1:(Nm+Ng),:);... % input and inhibition
                max(Wmg*Voltage(1:Nm, :), zeros(Nm, Ns))]; % excitation 
            
            % the ones not in refractory get incremented as per usual
            Voltage = Voltage + dt*(input-Voltage)./tau_x + [randn(Nm, Ns)*sqrt(dt); zeros(Ng, Ns)];
                        
            Voltage_history(:, :, count) = Voltage;
        end
        
        fprintf('------END------\n');

	fprintf('About to save\n');
        try
        	save([file_location, 'data.mat'], 'Voltage_history');
        catch
		continue;
	end
	fprintf('Successfully saved\n')
	end
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
    copyfile('Fopt_2MC_linear.m', file_location);
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
