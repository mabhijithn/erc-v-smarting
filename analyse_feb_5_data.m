clear;
close all;

if 0
% Location on your PC of the data (Feb-05)
main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-v-smarting-Feb-5-2020';

% Add the peak-detection function to the MATLAB path
% The toolbox is attached with the email
addpath('/home/abhijith/Documents/MATLAB/tools/pan_tompkin');

oddball_erc_file = 'triggerTest_05_02_2020.csv';
oddball_smarting_file = 'ercboard-v-smarting-trigger-test-1-[2020.02.05-15.15.50].gdf';

% Following files are with phones on airplance mode
%    oddball_erc_file = 'triggerTest2_05_02_2020';
%    oddball_smarting_file = 'ercboard-v-smarting-trigger-test-2-[2020.02.05-15.46.04].gdf';

EEG_smarting = pop_biosig(fullfile(main_fldr, oddball_smarting_file));
EEG_ercboard = csvread(fullfile(main_fldr, oddball_erc_file),2,0);
fs_smarting = EEG_smarting.srate;
fs_erc = 1000;

smarting_data = EEG_smarting.data(1,:); % Sine-wave at 10Hz
erc_sine_data = EEG_ercboard(:,3).*1000; % sine-wave recorded by ERC
erc_trigger_data = EEG_ercboard(:,2); % Audio triggers recorded by ERC


x_axis = 1:size(erc_trigger_data,1);
figure
plot(x_axis,erc_trigger_data(1:length(x_axis)));
title('Recorded Audio triggers: NOTE (Zoomin)');

% Extracting energy of the signal to find peaks
% Energy in a 60ms window extracted (the duration of a tone is 60ms)
win_size = 60e-3;
win_smpls = win_size*fs_erc;
strt = 1;
stp  = strt+win_smpls-1;
k = 1;
erc_energy_trigger = zeros(size(erc_trigger_data,1),1);
while(stp<size(erc_trigger_data,1))
    win_signal = erc_trigger_data(strt:stp);
    erc_energy_trigger(k) = sum(win_signal.^2)/win_smpls;
    strt = strt+1;
    stp = strt+win_smpls-1;
    k = k+1;
end



% Using an R-peak detector to pick the peaks 
% The function is attached in the email
plot_r_peak_stats = 0;
[qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(erc_energy_trigger,fs_erc,plot_r_peak_stats);
id = find(qrs_amp_raw>0.4);

% Index of oddball events in ERC data
sel_events = qrs_i_raw(id);
% Index of oddball events in Smarting data
smarting_events = find([EEG_smarting.event.type]==33025 | [EEG_smarting.event.type]==33026);

% Plotting the selected peaks. Overlayed on the energy of trigger waveform
% You can observe the peaks match.
x_axis = 1:size(erc_energy_trigger,1);
figure
plot(x_axis,erc_energy_trigger(1:length(x_axis)));
title('Energy in 60ms windows of the recorded trigger + Detected peaks');
hold on;
plot(sel_events, erc_energy_trigger(sel_events), '*');
end

