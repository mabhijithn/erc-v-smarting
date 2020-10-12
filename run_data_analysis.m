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
%     main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-new-bandpass-cutoff-July-31/data';
%     main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/1-Oct-2020';
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/11-Oct-2020';
    
    % ERC EEG file
    erc_file = 'alpha_trial_03.csv';
    erc_file = 'mariaAlpha03.csv';
    erc_file = 'MariaALpha_12_10_20__01.csv';
    
    
    % Smarting EEG file
    smarting_file = 'alpha-trial-2-1min-[2020.08.31-16.32.51].gdf';
    smarting_file = 'two-electrode-circuit-alpha-2-60s-[2020.10.01-16.08.56].gdf' %  Three-electrode recording; channel 2 from scalp
    smarting_file = 'alphat-test-1-[2020.10.12-16.00.03].gdf';
    
    % Read smarting data
    EEG_smarting = pop_biosig(fullfile(main_fldr, smarting_file));
    
    % Extract recorded channels
%     EEG_smarting = pop_select( EEG_smarting,'channel',{'Channel 1'});
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
    eeg_raw_struct.data = 1e3*erc_eeg_data';
    eeg_raw_struct.nbchan = 1;
    eeg_raw_struct.trials = 1;
    eeg_raw_struct.srate = fs_erc;
    eeg_raw_struct.pnts = size(eeg_raw_struct.data,2);
    eeg_raw_struct.times = round((1:eeg_raw_struct.pnts)*1e3/fs_erc);
    
    %     savepath = 'aug-31/figs';
    savepath = 'oct-11/figs';
    %% Calculate FFT and average  
    eeg_test = eeg_raw_struct;
    [avg_fft_open_1,avg_fft_open_2,avg_fft_close] = alpha_analyse(eeg_test, fs_new, dur);    
    ffts_erc = [avg_fft_open_1,avg_fft_close,avg_fft_open_2];
    
    eeg_test = EEG_smarting;
    [avg_fft_open_1,avg_fft_open_2,avg_fft_close] = alpha_analyse(eeg_test, fs_new, dur);
    ffts_smarting = [avg_fft_open_1,avg_fft_close,avg_fft_open_2];
    fftstoplot = struct;
    save_str = {'eye_open_1','eye_close_1','eye_open_2'};
    for j = 1:3
        for i = 1:2
            if i == 1
                ffts = ffts_erc;
                title_str = {'ERC Eyes Open #1 (60 sec)','ERC Eyes Closed (60 sec)','ERC Eyes Open #2 (60 sec)'};                
            else
                ffts = ffts_smarting;
                title_str = {'Smarting Eyes Open #1 (60 sec)','Smarting Eyes Closed (60 sec)','Smarting Eyes Open #2 (60 sec)'};                
            end
            fftstoplot(i).data = ffts(:,j);
            fftstoplot(i).fs = fs_new;
            fftstoplot(i).L = size(fftstoplot(i).data,1);
            fftstoplot(i).title = title_str{j};
            fftstoplot(i).savestr = save_str{j};
        end
        plot_fft_analysis(fftstoplot, 1, savepath);
    end
end

if 0    
    %% Time domain analysis of data
    
    smarting_file = 'alphat-test-1-[2020.10.12-16.00.03].gdf';
    erc_file = 'MariaALpha_12_10_20__01.csv';
    EEG_smarting = pop_biosig(fullfile(main_fldr, oddball_smarting_file));    
    EEG_smarting = pop_select( EEG_smarting,'channel',{'Channel 1'});
     %Index of audio trigger events in Smarting data
    smarting_events = find([EEG_smarting.event.type]==33027);
    
    EEG_smarting = pop_resample(EEG_smarting, fs_new);
    EEG_smarting = pop_eegfiltnew(EEG_smarting, [],0.5,[],true); %HP filtering (Filter DC)
    EEG_smarting = pop_eegfiltnew(EEG_smarting, [],20); %LP filtering
    
    for i = 1:numel(smarting_events)
        event_lats(i) = EEG_smarting.event(smarting_events(i)).latency;
    end
    
    EEG_ercboard = csvread(fullfile(main_fldr, oddball_erc_file),2,0);
    fs_erc = 1000;   
    
%     smarting_trigger_data = EEG_smarting.data(2,:); % EEG data
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
    for i = 1:numel(sel_events)
        eeg_raw_struct.event(i).type = 33026;
        eeg_raw_struct.event(i).latency = sel_events(i);
        eeg_raw_struct.event(i).urevent = i;
    end
    eeg_raw_struct.times = round((1:eeg_raw_struct.pnts)*1e3/fs_erc);
    eeg_raw_struct = pop_resample(eeg_raw_struct, fs_new);
    eeg_raw_struct = pop_eegfiltnew(eeg_raw_struct, [],0.5,[],true); %HP filtering (Filter DC)
    eeg_raw_struct = pop_eegfiltnew(eeg_raw_struct, [],20); %LP filtering
    %%
    % Add Smarting data as an extra channel in eeg_raw_struct
    eeg_raw_struct.nbchan = 2;
    if(eeg_raw_struct.pnts<EEG_smarting.pnts)
        minlen = eeg_raw_struct.pnts;
    else
        minlen = EEG_smarting.pnts;
        eeg_raw_struct.pnts = minlen;
        eeg_raw_struct.times = eeg_raw_struct.times(1:minlen);
    end
    eeg_raw_struct.data(2,1:minlen) = EEG_smarting.data(1,1:minlen);
    num_of_events = size(eeg_raw_struct.event,2);
    for i = 1:size(EEG_smarting.event,2)
        eeg_raw_struct.event(num_of_events+i) = EEG_smarting.event(i);
        eeg_raw_struct.event(num_of_events+i).urevent = num_of_events+i;
    end
    %%    
    eeg_raw_struct.data(1,:) = 1e3*eeg_raw_struct.data(1,:);
%     pop_eegplot(eeg_raw_struct);
    [c,lags] = xcorr(eeg_raw_struct.data(1,:),eeg_raw_struct.data(2,:));
    [mx,id] = max(c);
    lagval = lags(id);
    if(lagval<0)
        lagval = 0;
    end
    
    tmpdat = eeg_raw_struct.data(1,lagval+1:end);
    tmpdat = [tmpdat,zeros(1,lagval)];
    eeg_raw_struct.data(1,:) = tmpdat;
    
    dur = 5;
    smpls = dur*fs_new;
    strt = 1;
    stp = strt+smpls-1;
    x_time = (strt:stp)*(1/fs_new);
%     
    figure,
    while (stp<=size(eeg_raw_struct.data,2))
        plot(x_time,eeg_raw_struct.data(1,strt:stp));
        hold on;
        plot(x_time,eeg_raw_struct.data(2,strt:stp));
        xlabel('Time (s)');
        ylabel('Amplitude (uV)')
        hold off;
        legend('ERC','Smarting');
        axis([x_time(1) x_time(end) -70 70])
        strt = stp+1;
        stp = strt+smpls-1;
        x_time = (strt:stp)*(1/fs_new);
    end
    
    dur = 2;
    smpls = dur*fs_new;
    strt = 1;
    stp = strt+smpls-1;
    
    L = smpls;
      %calculate frequency bins with FFT
    df=fs_new/L; %frequency resolution
    sampleIndex = 0:L-1; %raw index for FFT plot
    f=sampleIndex*df; %x-axis index converted to frequencies 
   
    figure,
    while (stp+lagval<=size(eeg_raw_struct.data,2))
        str = sprintf('%3.2f to %3.2f sec',strt/fs_new,stp/fs_new);
        erc_sig = eeg_raw_struct.data(1,strt+lagval:stp+lagval);
        erc_sig_fft = fft(erc_sig);
        smarting_sig = eeg_raw_struct.data(2,strt:stp);     
        smarting_sig_fft = fft(smarting_sig);
        subplot(211)
        plot(f(1:floor(L/2)),abs(erc_sig_fft(1:floor(L/2))));
        ylabel('FFT Magnitude');
        xlabel('Frequency (Hz)');
        title(sprintf('ERC - %s',str));
        subplot(212)        
        plot(f(1:floor(L/2)),abs(smarting_sig_fft(1:floor(L/2))));
        ylabel('FFT Magnitude');
        xlabel('Frequency (Hz)');
        title(sprintf('Smarting- %s',str));
        strt = stp+1;
        stp = strt+smpls-1;
    end
end

if 0
    % Oddball experiment analysis
    
    % Location of data
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-new-bandpass-cutoff-July-31/data';
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/11-Oct-2020';
     
    % ERC EEG file
    oddball_erc_file = 'oddBall_01.csv';
    oddball_erc_file = 'MariaOddBall_12_10_20__01.csv';
    
    % Smarting EEG file
    oddball_smarting_file = 'ercboard-v-smarting-oddball-test-1-[2020.08.31-16.57.44].gdf';
    oddball_smarting_file = 'two-electrode-oddball-1-[2020.10.12-16.54.14].gdf';
    
    % Read ERC and Smarting data
    EEG_smarting = pop_biosig(fullfile(main_fldr, oddball_smarting_file));
    EEG_ercboard = csvread(fullfile(main_fldr, oddball_erc_file),2,0);
    fs_smarting = EEG_smarting.srate;
    fs_erc = 1000;
    
    fs_new = 120;
    
    smarting_data = EEG_smarting.data(1,:); % EEG data
    smarting_trigger_data = EEG_smarting.data(2,:); % EEG data
    erc_trigger_data = EEG_ercboard(:,3); % Audio triggers recorded by ERC
    
    peak_thr = 0.4; % For trial-1
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