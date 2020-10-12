clc;
clear;
close all;
if 1
    % Point-by-point analysis
    % Compare EEG recorded by Smarting and ERC point-by-point
    eeglab;
    fs_new = 120;
%     main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-board-July-28-2020/data';        
%     oddball_erc_file = 'ERC_alpha-trial.csv';
%     oddball_smarting_file = 'alpha-trial-1-[2020.07.28-16.22.34].gdf';    
%     oddball_erc_file = 'ERC_alpha-trial2.csv';
%     oddball_smarting_file = 'alpha-trial-2-[2020.07.28-16.44.11].gdf';
    
        %% August-31 data
%     main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-new-bandpass-cutoff-July-31/data';
%     oddball_erc_file = 'alpha_trial_03.csv';
%     oddball_smarting_file = 'alpha-trial-2-1min-[2020.08.31-16.32.51].gdf';
    
    %% October 1 data
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/1-Oct-2020';
    oddball_erc_file = 'mariaAlpha03.csv';
    oddball_smarting_file = 'two-electrode-circuit-alpha-2-60s-[2020.10.01-16.08.56].gdf' %  Three-electrode recording; channel 2 from scalp
    
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/11-Oct-2020';
    oddball_erc_file = 'MariaALpha_12_10_20__01.csv';
    oddball_smarting_file = 'alphat-test-1-[2020.10.12-16.00.03].gdf';
    %%
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
%     for i = 1:numel(smarting_events)-1
%         eeg_smarting = EEG_smarting.data(1,event_lats(i):event_lats(i+1));
%         len_smarting = numel(eeg_smarting);
%         eeg_erc = eeg_raw_struct.data(1,sel_events(i):2:sel_events(i+1),1);
%         len_erc = numel(eeg_erc);
%         if(len_erc<len_smarting)
%             minlen = len_erc;
%         else
%             minlen = len_smarting;
%         end
%         figure
%         plot(1000.*eeg_erc(1:minlen))
%         hold on
%         plot(eeg_smarting(1:minlen));
%         hold off;
%     end
end

if 0
    % Alpha wave analysis
    eeglab;
    fs_new = 120;
    %% July-28 data
%     main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-board-July-28-2020/data';
    % Following files are with phones on airplance mode
%     oddball_erc_file = 'ERC_eyeblink-trial.csv';
%     oddball_smarting_file = 'eyeblink-trial-[2020.07.28-15.42.29].gdf';     
%     oddball_erc_file = 'ERC_alpha-trial.csv';
%     oddball_smarting_file = 'alpha-trial-1-[2020.07.28-16.22.34].gdf';
%     
%     oddball_erc_file = 'ERC_alpha-trial2.csv';
%     oddball_smarting_file = 'alpha-trial-2-[2020.07.28-16.44.11].gdf';
    %% August-31 data
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-new-bandpass-cutoff-July-31/data';
    oddball_erc_file = 'alpha_trial_03.csv';
    oddball_smarting_file = 'alpha-trial-2-1min-[2020.08.31-16.32.51].gdf';
    %%
    EEG_smarting = pop_biosig(fullfile(main_fldr, oddball_smarting_file));
    
    EEG_smarting = pop_select( EEG_smarting,'channel',{'Channel 1'});
%     EEG_smarting = pop_resample(EEG_smarting, fs_new);
%     EEG_smarting = pop_eegfiltnew(EEG_smarting, [],0.5,[],true); %HP filtering
%     EEG_smarting = pop_eegfiltnew(EEG_smarting, [],20); %LP filtering
%     smarting_eeg_data_filtered = EEG_smarting.data(1,:)'; % EEG data
    
    EEG_ercboard = csvread(fullfile(main_fldr, oddball_erc_file),2,0);
    fs_smarting = EEG_smarting.srate;
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
    plot_fft_analysis(fftstoplot, 1, savepath);
                
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
    plot_fft_analysis(fftstoplot, 1, savepath);

end

if 0
    % Oddball experiment analysis
    
%     main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-board-July-28-2020/data';
    
    % Add the peak-detection function to the MATLAB path
    % The toolbox is attached with the email
    addpath('/home/abhijith/Documents/MATLAB/toolboxes/pan_tompkin');
    
%     oddball_smarting_file = 'ercboard-v-smarting-oddball-test-1-[2020.07.28-17.00.19].gdf';
%     oddball_erc_file = 'ERC_oddBall_trial.csv';
    
%     oddball_smarting_file = 'ercboard-v-smarting-oddball-test-2-without-sound-noise-[2020.07.28-17.34.05].gdf';
%     oddball_erc_file = 'ERC_oddBall_trial2.csv';
    
    %%
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-new-bandpass-cutoff-July-31/data';
    oddball_erc_file = 'oddBall_01.csv';
    oddball_smarting_file = 'ercboard-v-smarting-oddball-test-1-[2020.08.31-16.57.44].gdf';
    
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
    
    % Analyse Smarting data
    oddball_analyse(EEG_smarting, 1, 'Smarting Data');
end

if 0
    % Oddball experiment analysis
    eeglab;
    close all;
%     main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/erc-board-July-28-2020/data';
    
    % Add the peak-detection function to the MATLAB path
    % The toolbox is attached with the email
    addpath('/home/abhijith/Documents/MATLAB/toolboxes/pan_tompkin');
    
%     oddball_smarting_file = 'ercboard-v-smarting-oddball-test-1-[2020.07.28-17.00.19].gdf';
%     oddball_erc_file = 'ERC_oddBall_trial.csv';
    
%     oddball_smarting_file = 'ercboard-v-smarting-oddball-test-2-without-sound-noise-[2020.07.28-17.34.05].gdf';
%     oddball_erc_file = 'ERC_oddBall_trial2.csv';
    
    
    EEG_smarting = pop_biosig(fullfile(main_fldr, oddball_smarting_file));
    EEG_smarting = pop_select( EEG_smarting,'channel',{'Channel 1'});
    
    EEG_ercboard = csvread(fullfile(main_fldr, oddball_erc_file),2,0);
    erc_trigger_data = EEG_ercboard(:,3); % Audio triggers recorded by ERC
    fs_smarting = EEG_smarting.srate;
    fs_erc = 1000;
    fs_new = 120;
    
    smarting_data = EEG_smarting.data(1,:); % EEG data
    
    peak_thr = 0.315; % For trial-2
    sel_events = extract_trigger_events(erc_trigger_data, fs_erc, peak_thr);
    
    %% Prepare EEGLAB structure from ERC data
    EEG_erc = EEG_smarting;
    EEG_erc.nbchan = 2;
    EEG_erc.trials = 1;
    EEG_erc.data = [EEG_ercboard(:,2)';zeros(1,size(EEG_ercboard,1))];
    EEG_erc.srate = fs_erc;    
    EEG_erc.pnts = size(EEG_erc.data,2);
    EEG_erc.times = (0:EEG_erc.pnts-1).*(1/EEG_erc.srate);
    
    % Index of oddball events in Smarting data
    smarting_events = find([EEG_smarting.event.type]==33025 | [EEG_smarting.event.type]==33026); 
    % The latencies of events are replaced by those extracted from audio
    % triggers
    for i = 1:numel(sel_events)
        if(EEG_smarting.event(smarting_events(i)).type==33025)
            EEG_erc.event(i).type = 33027;
        else
            EEG_erc.event(i).type = 33028;
        end
        EEG_erc.event(i).latency = sel_events(i);
        EEG_erc.event(i).urevent = i;
    end
    
    %% Add smarting data as an extra channel
    EEG_erc = pop_resample(EEG_erc, fs_smarting);
    
    if(EEG_erc.pnts<EEG_smarting.pnts)
        minlen = EEG_erc.pnts;
    else
        minlen = EEG_smarting.pnts;
        EEG_erc.pnts = minlen;
        EEG_erc.times = EEG_erc.times(1:minlen);
    end
    EEG_erc.data(2,1:minlen) = EEG_smarting.data(1,1:minlen);
    num_of_events = size(EEG_erc.event,2);
    for i = 1:size(EEG_smarting.event,2)
        EEG_erc.event(num_of_events+i) = EEG_smarting.event(i);
        EEG_erc.event(num_of_events+i).urevent = num_of_events+i;
    end
    % Scaling for plotting    
    EEG_erc.data(1,:) = double(1e3*EEG_erc.data(1,:));
    
    pop_eegplot(EEG_erc);
 end
