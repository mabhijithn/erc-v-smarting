clc;
clear;
close all;
% Single-trial classification of oddball events using DDTBOX

% Location on your PC of the data (Feb-12)
main_fldr = 'C:\Users\abhijith\Documents\MEGA\MEGAsync\PhD\2019\masters_thesis_micas\ERC_vs_Smarting\Feb-12';

% Add the peak-detection function to the MATLAB path
% The toolbox is attached with the email
addpath('C:\Users\abhijith\Documents\MATLAB\tools\pan_tompkin');

 % Following files are with phones on airplance mode
oddball_erc_file = 'oddballTest1_2020_02_15.csv';
oddball_smarting_file = 'ercboard-v-smarting-oddball-test-1-[2020.02.12-15.21.46].gdf';


EEG_smarting = pop_biosig(fullfile(main_fldr, oddball_smarting_file));

%% Convertion of EEGLAB epoched data to DDTBOX format for classification
save_directory  = 'C:\Users\abhijith\Documents\MATLAB\erc_v_smarting';
save_filename = 'Feb12_smarting_oddball_DDTBOX_frmt.mat';
EEGB = pop_epoch( EEG_smarting, {  '33025' ,'33026' }, [-0.2  1.0], 'newname', 'GDF file epochs', 'epochinfo', 'yes'); %025 Extract Background
EEGD = pop_epoch( EEG_smarting, {  '33026'  }, [-0.2  1.0], 'newname', 'GDF file epochs', 'epochinfo', 'yes'); %026 Extract Target High tone

events_by_cond{1, 1} = [33025];%, 33024];
events_by_cond{1, 2} = [33026];%, 33024];

dd_convert_eeg_data(EEGB, events_by_cond, save_directory, save_filename);

