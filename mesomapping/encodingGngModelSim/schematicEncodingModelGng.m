%This encoding toy model code is modified from: /Users/jp3025/Library/CloudStorage/GoogleDrive-jp3025@princeton.edu/My Drive/Mover_shared/codes/Figure1_SchematicModel.m

%% Parameters
N = 120;  % Total neurons
t_steps = 2500;  % Time steps per trial
num_sessions = 20;  % Total learning sessions
num_trials_per_session = 20;  % 10 go trials and 10 no-go trials per session
learning_rate = linspace(0.1, 1, num_sessions);  % Scaling from early to late learning

% Parameters for promiscuous responses
initial_overlap_scale = 0.8;  % Higher starting overlap for more promiscuity
final_overlap_scale = 0.3;    % Set a minimum overlap to retain promiscuity
response_jitter = 0.3;        % Small jitter to allow for overlap without large trial-to-trial variance

% Base response values for go and no-go tones
base_go_response = 5;  
base_no_go_response = 5;

% Generate variable lags for Cue and Action responses
t_lags_all_cue = poissrnd((1:N) * 2);  % Poisson-distributed lags for cue response
tmp = sort(t_lags_all_cue);
t_lags_cue = [tmp(1:2:end) tmp(2:2:end)] + 500;  % Add base delay

t_lags_all_action = poissrnd((1:N) * 9);  % Larger range for action response
tmp = sort(t_lags_all_action);
t_lags_action = [tmp(1:2:end) tmp(2:2:end)] + 700;  % Add base delay

% figure; hold on; 
% plot(t_lags_cue); 
% plot(t_lags_action); 
% grid on
% title("Time lag for cue and action")
% xlabel("Neuron")
% ylabel("Time(ms)")
% hold off;
% print(fullfile('/Volumes/buschman/Rodent Data/dualImaging_parkj/collectFigure', 'timeLag_cue_action'), '-bestfit', '-vector', '-dpdf'); 

% Create separate kernels for cue and action responses for each neuron
kernels = zeros(N * 2, 2e3);  % Double size to store separate kernels for cue and action

for ii = 1:N
    % Cue kernel
    [kernel_rise] = TNC_CreateGaussian(1e3, 25 + rand(1) * 40, 2e3, 1);
    kernel_rise = kernel_rise / max(kernel_rise);
    [kernel] = TNC_CreateGaussian(1e3, 50 + rand(1) * 200, 2e3, 1);
    kernel = kernel / max(kernel);
    kernel(1:1e3) = kernel_rise(1:1e3);
    kernels(ii, :) = kernel;  % Assign to cue kernel set

    % Action kernel
    [kernel_rise] = TNC_CreateGaussian(1e3, rand(1) * ((t_lags_action(ii) - 600) / 4), 2e3, 1);
    kernel_rise = kernel_rise / max(kernel_rise);
    [kernel] = TNC_CreateGaussian(1e3, 50 + rand(1) * 300, 2e3, 1);
    kernel = kernel / max(kernel);
    kernel(1:1e3) = kernel_rise(1:1e3);
    kernels(ii + N, :) = kernel;  % Assign to action kernel set
end

% Set an equal ratio of go and no-go neurons
go_tone_neurons = N / 2; 
no_go_tone_neurons = N - go_tone_neurons;

% Randomly select 50% of neurons to be action-responsive
action_tuned_neurons = randperm(N, round(N * 0.5));

% Initialize selector.psth as a cell array
selector.psth = cell(1, num_sessions);  % Cell array to hold data for each session

%% Main Loop over Sessions
for session = 1:num_sessions
    % Update overlap and amplitude based on learning progression
    overlap_scale = max(1 - learning_rate(session) * (initial_overlap_scale - final_overlap_scale), final_overlap_scale);
    amplitude_scale = 1 + (learning_rate(session) * 0.7);  % Scales cue response over learning
    action_selectivity = learning_rate(session);  % Higher selectivity in later sessions

    % Initialize the 3D matrix for the current session (120 neurons x 2500 time steps x 20 trials)
    session_data = zeros(N, t_steps, num_trials_per_session);

    % Loop over trials within each session (10 go, 10 no-go)
    for trial = 1:num_trials_per_session
        is_go_trial = mod(trial, 2) == 1;  % Odd trials as Go, even as No-Go

        for ii = 1:N
            % Cue Response (Unchanged)
            if ii <= go_tone_neurons
                go_tone_response = (base_go_response + rand(1) * response_jitter) * amplitude_scale;
                no_go_tone_response = (base_no_go_response + rand(1) * response_jitter) * overlap_scale;
            else
                go_tone_response = (base_go_response + rand(1) * response_jitter) * overlap_scale;
                no_go_tone_response = (base_no_go_response + rand(1) * response_jitter) * amplitude_scale;
            end
            
            % Initialize delta for cue response and assign based on trial type
            delta = zeros(1, t_steps);
            if is_go_trial
                delta(t_lags_cue(ii)) = go_tone_response;
            else
                delta(t_lags_cue(ii)) = no_go_tone_response;
            end

            % Convolve with cue-specific kernel and store in the session data
            session_data(ii, :, trial) = conv(delta, kernels(ii, :), 'same');

            % Action Response (Modified to Gain Selectivity Across Sessions)
            if ismember(ii, action_tuned_neurons)
                % Fixed magnitude for action response across all sessions
                fixed_action_response = base_go_response + rand(1) * response_jitter;

                % Initialize delta for action response
                delta = zeros(1, t_steps);
                if is_go_trial || rand < (1 - action_selectivity)
                    % Start with promiscuous response; gradually focus on go trials
                    delta(t_lags_action(ii)) = fixed_action_response;
                elseif ~is_go_trial && rand < (1 - action_selectivity)
                    % Reduce response in no-go trials over sessions
                    delta(t_lags_action(ii)) = fixed_action_response * (1 - action_selectivity);
                end

                % Convolve with action-specific kernel and add to the session data
                session_data(ii, :, trial) = session_data(ii, :, trial) + conv(delta, kernels(ii + N, :), 'same');
            end
        end
    end

    % Store the session data in the selector.psth cell array
    selector.psth{session} = session_data;
    fprintf("Completed computing the session #%d\n", session)
end

%figure; imagesc(selector.psth{end}(:, :, 20));
save(fullfile('/Users/jp3025/Documents/codes/MatlabNeuralDataPipeline/mesomapping/encodingGngModelSim/simData', 'simData'), 'selector', '-mat'); 
