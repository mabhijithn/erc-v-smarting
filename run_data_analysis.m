clc;
clear;
close all;

% Add the peak-detection function to the MATLAB path
addpath('/home/abhijith/Documents/MATLAB/toolboxes/pan_tompkin');

if 1
    % Alpha analysis
    eeglab; close all;
    
    % Duration of eyes-closed in the experiment
    dur = 60; % in seconds
    
    % Downsample frequency
    fs_new = 120;
    
    % Location of data
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-new-bandpass-cutoff-July-31/data';
    
    % ERC EEG file
    erc_file = 'alpha_trial_03.csv';
    
    % Smarting EEG file
    smarting_file = 'alpha-trial-2-1min-[2020.08.31-16.32.51].gdf';
    
    % Read smarting data
    EEG_smarting = pop_biosig(fullfile(main_fldr, smarting_file));
    
    % Extract recorded channels
    EEG_smarting = pop_select( EEG_smarting,'channel',{'Channel 1'});
    fs_smarting = EEG_smarting.srate;
    
    % Read ERC data
    EEG_ercboard = csvread(fullfile(main_fldr, erc_file),2,0);
    fs_erc = 1000;
    
    erc_trigger_data = EEG_ercboard(:,3); % Audio triggers recorded by ERC
    erc_eeg_data = EEG_ercboard(:,2); %EEG Data
    
    peak_thr = 0.3;
    sel_events = extract_trigger_events(erc_trigger_data, fs_erc, peak_thr);
    
    %% Create an EEGLAB structure for ERC data
    eeg_raw_struct = EEG;
    eeg_raw_struct.data = erc_eeg_data';
    eeg_raw_struct.nbchan = 1;
    eeg_raw_struct.trials = 1;
    eeg_raw_struct.srate = fs_erc;
    eeg_raw_struct.pnts = size(eeg_raw_struct.data,2);
    eeg_raw_struct.times = round((1:eeg_raw_struct.pnts)*1e3/fs_erc);
    
    %% Calculate FFT and average  
    eeg_test = eeg_raw_struct;
    [avg_fft_open_1,avg_fft_open_2,avg_fft_close] = alpha_analyse(eeg_test, fs_new, dur);
    
    savepath = 'aug-31/figs';
    ffts = [avg_fft_open_1,avg_fft_close,avg_fft_open_2];
    title_str = {'ERC Eyes Open #1 (60 sec)','ERC Eyes Closed (60 sec)','ERC Eyes Open #2 (60 sec)'};
    save_str = {'erc_eye_open_1','erc_eye_close_1','erc_eye_open_2'};
    for i = 1:3
        fftstoplot(i).data = ffts(:,i);
        fftstoplot(i).fs = fs_new;
        fftstoplot(i).L = size(fftstoplot(i).data,1);
        fftstoplot(i).title = title_str{i};
        fftstoplot(i).savestr = save_str{i};
    end
    plot_fft_analysis(fftstoplot, 0, savepath);
                
    eeg_test = EEG_smarting;
    [avg_fft_open_1,avg_fft_open_2,avg_fft_close] = alpha_analyse(eeg_test, fs_new, dur);
    ffts = [avg_fft_open_1,avg_fft_close,avg_fft_open_2];
    title_str = {'Smarting Eyes Open #1 (60 sec)','Smarting Eyes Closed (60 sec)','Smarting Eyes Open #2 (60 sec)'};
    save_str = {'smarting_eye_open_1','smarting_eye_close_1','smarting_eye_open_2'};
    for i = 1:3
        fftstoplot(i).data = ffts(:,i);
        fftstoplot(i).fs = fs_new;
        fftstoplot(i).L = size(fftstoplot(i).data,1);
        fftstoplot(i).title = title_str{i};
        fftstoplot(i).savestr = save_str{i};
    end
    plot_fft_analysis(fftstoplot, 0, savepath);
end

if 1
    % Oddball experiment analysis
    
    % Location of data
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-new-bandpass-cutoff-July-31/data';
    
    % ERC EEG file
    oddball_erc_file = 'oddBall_01.csv';
    
    % Smarting EEG file
    oddball_smarting_file = 'ercboard-v-smarting-oddball-test-1-[2020.08.31-16.57.44].gdf';
    
    % Read ERC and Smarting data
    EEG_smarting = pop_biosig(fullfile(main_fldr, oddball_smarting_file));
    EEG_ercboard = csvread(fullfile(main_fldr, oddball_erc_file),2,0);
    fs_smarting = EEG_smarting.srate;
    fs_erc = 1000;
    
    fs_new = 120;
    
    smarting_data = EEG_smarting.data(1,:); % EEG data
    smarting_trigger_data = EEG_smarting.data(2,:); % EEG data
    erc_trigger_data = EEG_ercboard(:,3); % Audio triggers recorded by ERC
    
    peak_thr = 0.25; % For trial-1
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
    
    % Analyse ERC data
    oddball_analyse(EEG_erc, 1, 'ERC Data');
    
    % Analyse Smarting data
    oddball_analyse(EEG_smarting, 1, 'Smarting Data');
end