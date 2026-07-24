%% Sign-flipped permutation
%% Paths
projectFolder = [];
rawDataFolder = [projectFolder,filesep,'Data'];
addpath(genpath('[]/eeglab2026.0.0'));
addpath([]); % matplotlib, load MAGMA as the color palette
figureFolder = [];
%% Load the data
dataFolder = [projectFolder,filesep,'Data',filesep,'unfold'];
if ERP == 1
    load([dataFolder,filesep,'correlationERP.mat'])
elseif POWER == 1
    load([dataFolder,filesep,'correlationTF.mat'])
end
%%
rsa_corr = similarity;
[nSubjects, nElectrodes] = size(rsa_corr);
nPerm = 10000;
%% Fisher transform
% Convert r values to z values
rsa_z = atanh(rsa_corr);
%% Initialize outputs
mean_z = zeros(1,nElectrodes);
mean_r = zeros(1,nElectrodes);
p_values = zeros(1,nElectrodes);
%% Permutation test for each electrode
for elec = 1:nElectrodes
    % Values across subjects for this electrode
    data = rsa_z(:,elec);
    % Observed statistic
    obs_mean = mean(data);
    mean_z(elec) = obs_mean;
    % Null distribution
    perm_means = zeros(nPerm,1);
    for p = 1:nPerm
        % Randomly flip sign for each subject
        signs = randi([0 1], nSubjects, 1);
        signs(signs==0) = -1;
        % Apply sign flipping
        shuffled_data = data .* signs;
        % Store mean
        perm_means(p) = mean(shuffled_data);
    end
    % One-sided p-value:
    % H0: mean correlation <= 0
    p_values(elec) = ...
        (sum(perm_means >= obs_mean)+1)/(nPerm+1);
end
%% Convert mean Fisher-z back to r for visualization
mean_r = tanh(mean_z);
%% Display results
results = table((1:nElectrodes)', ...
                mean_r', ...
                mean_z', ...
                p_values', ...
                'VariableNames',...
                {'Electrode','Mean_r','Mean_Fisher_z','p_value'});

disp(results)
%% Optional: significance threshold
alpha = 0.05;
sig_electrodes = p_values < alpha;
disp('Significant electrodes:')
disp(find(sig_electrodes))