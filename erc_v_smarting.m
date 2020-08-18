clc;
clear;
close all;
main_fldr = 'C:\Users\abhijith\Downloads\Data\erc-aquisition-v-smarting';
%%
erc_board = 'exp2_synthetic.csv';
smarting_data_name = 'synthetic-signals-test-2-[2019.12.11-14.25.34].gdf';

EEG_smarting = pop_biosig(fullfile(main_fldr, smarting_data_name));
EEG_ercboard = csvread(fullfile(main_fldr, erc_board),2,0);

smarting_data = EEG_smarting.data(1,:);
erc_data = EEG_ercboard(:,2);
fs_smarting = 500;
fs_erc = 1000;

figure
plot(smarting_data);
title('Smarting');

figure
plot(erc_data);
title('ERC');

smarting_start = 22527;
erc_start = 6149;

figure
plot(erc_data(erc_start:2:end));
hold on;
plot(smarting_data(smarting_start:end));
hold off;
legend('ERC','Smarting');

%% Eyeblink Test
eyeblink_smarting_file = 'eyeblink-test-1-[2019.12.11-15.08.35].gdf';
eyeblink_erc_file = 'eyeBlinkTest1.csv';

EEG_smarting = pop_biosig(fullfile(main_fldr, eyeblink_smarting_file));
EEG_ercboard = csvread(fullfile(main_fldr, eyeblink_erc_file),2,0);

smarting_data = EEG_smarting.data(1,:);
erc_data = EEG_ercboard(:,2);
fs_smarting = 500;
fs_erc = 1000;

figure
plot(erc_data(erc_start:2:end));
hold on;
plot(smarting_data(smarting_start:end));
hold off;
legend('ERC','Smarting');

%% Alpha Test
alpha_smarting_file = 'alpha-test-1-[2019.12.11-15.12.57].gdf';
alpha_erc_file = 'alphaTesst1.csv';
