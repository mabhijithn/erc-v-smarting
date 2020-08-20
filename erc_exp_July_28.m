clc;
clear;
close all;
if 1
    % Alpha wave analysis
    eeglab;
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-board-July-28-2020/data';
    % Following files are with phones on airplance mode
%     oddball_erc_file = 'ERC_eyeblink-trial.csv';
%     oddball_smarting_file = 'eyeblink-trial-[2020.07.28-15.42.29].gdf';
    fs_new = 120;
    
    oddball_erc_file = 'ERC_alpha-trial.csv';
    oddball_smarting_file = 'alpha-trial-1-[2020.07.28-16.22.34].gdf';
    
    oddball_erc_file = 'ERC_alpha-trial2.csv';
    oddball_smarting_file = 'alpha-trial-2-[2020.07.28-16.44.11].gdf';
    
    EEG_smarting = pop_biosig(fullfile(main_fldr, oddball_smarting_file));
    
    EEG_smarting = pop_select( EEG_smarting,'channel',{'Channel 1'});
    EEG_smarting = pop_resample(EEG_smarting, fs_new);
    EEG_smarting = pop_eegfiltnew(EEG_smarting, [],0.5,[],true); %HP filtering
    EEG_smarting = pop_eegfiltnew(EEG_smarting, [],20); %LP filtering
    smarting_eeg_data_filtered = EEG_smarting.data(1,:)'; % EEG data
    
    EEG_ercboard = csvread(fullfile(main_fldr, oddball_erc_file),2,0);
    fs_smarting = EEG_smarting.srate;
    fs_erc = 1000;   
    
%     smarting_trigger_data = EEG_smarting.data(2,:); % EEG data
    erc_trigger_data = EEG_ercboard(:,3); % Audio triggers recorded by ERC
    erc_eeg_data = EEG_ercboard(:,2); %EEG Data
    
    %% Create an EEGLAB structure for ERC data
    eeg_raw_struct = EEG;
    eeg_raw_struct.data = erc_eeg_data';
    eeg_raw_struct.nbchan = 1;
    eeg_raw_struct.trials = 1;
    eeg_raw_struct.srate = fs_erc;
    eeg_raw_struct.pnts = size(eeg_raw_struct.data,2);
    eeg_raw_struct.times = round((1:eeg_raw_struct.pnts)*1e3/fs_erc);
%% Create Band-pass filter
%     Fs = 120; %1000;
%     fl = 20;
%     fh = 0.5;
%     [ BP_equirip ] = cnstr_bpfilter(Fs, fl, fh);
%% Filter EEG data    
%     erc_eeg_resampled = resample(erc_eeg_data,120,fs_erc);
%     erc_eeg_data_filtered = filtfilt(BP_equirip.numerator, 1, double(erc_eeg_resampled));
    
%     smarting_eeg_resampled = resample(double(smarting_data'),120,fs_smarting);
%     smarting_eeg_data_filtered = filtfilt(BP_equirip.numerator, 1, double(smarting_eeg_resampled));
%% Calculate FFT and average        
    
    [avg_fft_open_1,avg_fft_open_2,avg_fft_close] = alpha_analyse(eeg_raw_struct, fs_new);
    
    %calculate frequency bins with FFT
    df=fs_new/L; %frequency resolution
    sampleIndex = 0:L-1; %raw index for FFT plot
    f=sampleIndex*df; %x-axis index converted to frequencies
    plot(f(1:floor(L/2)),abs(avg_fft_open_1(1:floor(L/2))));
    ylabel('FFT Magnitude (Avg)');
    xlabel('Frequency (Hz)');
    title('Eyes Open #1 (60 sec)');
%     saveas(gcf, 'july-28/fig/eyes_open_erc_60s_1.jpg');
    
    figure,
    plot(f(1:L/2),abs(avg_fft_close(1:L/2)));
    ylabel('FFT Magnitude (Avg)');
    xlabel('Frequency (Hz)');
    title('Eyes Closed (60 sec)');
%     saveas(gcf, 'july-28/fig/eyes_closed_erc_60s_1.jpg');
    
    figure
    plot(f(1:L/2),abs(avg_fft_open_2(1:L/2)));
    ylabel('FFT Magnitude (Avg)');
    xlabel('Frequency (Hz)');
    title('Eyes Open #2 (60 sec)');
%     saveas(gcf, 'july-28/fig/eyes_open_erc_60s_2.jpg');
end

if 0
    % Oddball experiment analysis
    
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-board-July-28-2020/data';
    
    % Add the peak-detection function to the MATLAB path
    % The toolbox is attached with the email
    addpath('/home/abhijith/Documents/MATLAB/toolboxes/pan_tompkin');
    
%     oddball_smarting_file = 'ercboard-v-smarting-oddball-test-1-[2020.07.28-17.00.19].gdf';
%     oddball_erc_file = 'ERC_oddBall_trial.csv';
    
    oddball_smarting_file = 'ercboard-v-smarting-oddball-test-2-without-sound-noise-[2020.07.28-17.34.05].gdf';
    oddball_erc_file = 'ERC_oddBall_trial2.csv';
    
    EEG_smarting = pop_biosig(fullfile(main_fldr, oddball_smarting_file));
    EEG_ercboard = csvread(fullfile(main_fldr, oddball_erc_file),2,0);
    fs_smarting = EEG_smarting.srate;
    fs_erc = 1000;
    fs_new = 120;
    
    smarting_data = EEG_smarting.data(1,:); % EEG data
    smarting_trigger_data = EEG_smarting.data(2,:); % EEG data
    erc_trigger_data = EEG_ercboard(:,3); % Audio triggers recorded by ERC
    x_axis = 1:size(erc_trigger_data,1);
    figure
    plot(x_axis,erc_trigger_data(1:length(x_axis)));
    title('Recorded Audio triggers: NOTE (Zoomin)');

    peak_thr = 0.3; % For trial-1
%     peak_thr = 0.315; % For trial-2
    sel_events = extract_trigger_events(erc_trigger_data, fs_erc, peak_thr);

    % Analyse Smarting data
    EEG_smarting.data = -1*EEG_smarting.data;
    oddball_analyse(EEG_smarting, 1, 'Smarting Data'); 
    
    % Index of oddball events in Smarting data
    smarting_events = find([EEG_smarting.event.type]==33025 | [EEG_smarting.event.type]==33026); 
    
    % Prepare EEGLAB structure from ERC data
    EEG_erc = EEG_smarting;
    EEG_erc.data = EEG_ercboard(:,2)';
    EEG_erc.srate = fs_erc;
    EEG_erc.event = EEG_smarting.event;
    EEG_erc.pnts = size(EEG_erc.data,1);
    EEG_erc.times = (0:EEG_erc.pnts-1).*(1/EEG_erc.srate);
    
    % The latencies of events are replaced by those extracted from audio
    % triggers
    for i = 1:length(smarting_events)
        EEG_erc.event(smarting_events(i)).latency =  sel_events(i);
    end    
%    
%     EEG_erc.event(smarting_events(1)).latency = sel_events(1);
%     for i = 2:length(smarting_events)
%         EEG_erc.event(smarting_events(i)).latency =  EEG_erc.event(smarting_events(i-1)).latency+1000;
%     end   
    
%     strt = EEG_erc.event(smarting_events(1)).latency;
%     stp = EEG_erc.event(smarting_events(2)).latency;
%     figure,
%     subplot(211)
%     plot(EEG_erc.data(1,strt:stp-1));
%     strt = EEG_smarting.event(smarting_events(1)).latency;
%     stp = EEG_smarting.event(smarting_events(2)).latency; 
%     subplot(212)
%     plot(EEG_smarting.data(1,strt:stp-1));
%     
%     for i = 1:2:length(smarting_events)-2
% %         figure
%         strt = EEG_erc.event(smarting_events(i)).latency;
%         stp = EEG_erc.event(smarting_events(i+2)).latency;
% %         plot(1000.*EEG_erc.data(1,strt:10:stp-1))
%         erc_win_data = 1000.*EEG_erc.data(1,strt:10:stp-1);
% %         hold on;
%         strt = EEG_smarting.event(smarting_events(i)).latency;
%         stp = EEG_smarting.event(smarting_events(i+2)).latency; 
% %         plot(EEG_smarting.data(1,strt:stp-1));
%         smarting_win_data = EEG_smarting.data(1,strt:stp-1);
%         corr_sig = xcorr(erc_win_data, smarting_win_data(1:length(erc_win_data)));
%         plot(corr_sig);
%     end    
%     
    % Analyse ERC data
    oddball_analyse(EEG_erc, 1, 'ERC Data');
end
