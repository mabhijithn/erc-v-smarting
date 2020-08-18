if 1
    % Location on your PC of the data (Feb-05)
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/Feb-12';

    % Add the peak-detection function to the MATLAB path
    % The toolbox is attached with the email
    addpath('/home/abhijith/Documents/MATLAB/toolboxes/pan_tompkin');
    
    % Following files are with phones on airplance mode
    oddball_erc_file = 'oddballTest1_2020_02_15.csv';
    oddball_smarting_file = 'ercboard-v-smarting-oddball-test-1-[2020.02.12-15.21.46].gdf';
   
   
    EEG_smarting = pop_biosig(fullfile(main_fldr, oddball_smarting_file));
    EEG_ercboard = csvread(fullfile(main_fldr, oddball_erc_file),2,0);
    fs_smarting = EEG_smarting.srate;
    fs_erc = 5000;
    fs_new = 500;
    
    smarting_data = EEG_smarting.data(1,:); % EEG data
    smarting_trigger_data = EEG_smarting.data(2,:); % EEG data
    erc_trigger_data = EEG_ercboard(:,3); % Audio triggers recorded by ERC   
   
    x_axis = 1:size(erc_trigger_data,1);
    figure
    plot(x_axis,erc_trigger_data(1:length(x_axis)));
    title('Recorded Audio triggers: NOTE (Zoomin)');

    peak_thr = 0.4;
    sel_events = extract_trigger_events (erc_trigger_data, fs_erc, peak_thr);
    
    % Plotting the selected peaks. Overlayed on the energy of trigger waveform
    % You can observe the peaks match.
%     x_axis = 1:size(erc_energy_trigger,1);
%     figure
%     plot(x_axis,erc_energy_trigger(1:length(x_axis)));
%     title('Energy in 60ms windows of the recorded trigger + Detected peaks');
%     hold on;
%     plot(sel_events, erc_energy_trigger(sel_events), '*');
    
    % Analyse Smarting data
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
    
    strt = EEG_erc.event(smarting_events(1)).latency;
    stp = EEG_erc.event(smarting_events(2)).latency;
    figure,
    subplot(211)
    plot(EEG_erc.data(1,strt:stp-1));
    strt = EEG_smarting.event(smarting_events(1)).latency;
    stp = EEG_smarting.event(smarting_events(2)).latency; 
    subplot(212)
    plot(EEG_smarting.data(1,strt:stp-1));
    
    for i = 1:2:length(smarting_events)-1
%         figure
        strt = EEG_erc.event(smarting_events(1)).latency;
        stp = EEG_erc.event(smarting_events(2)).latency;
%         plot(1000.*EEG_erc.data(1,strt:10:stp-1))
        erc_win_data = 1000.*EEG_erc.data(1,strt:10:stp-1);
        hold on;
        strt = EEG_smarting.event(smarting_events(1)).latency;
        stp = EEG_smarting.event(smarting_events(2)).latency; 
%         plot(EEG_smarting.data(1,strt:stp-1));
        smarting_win_data = EEG_smarting.data(1,strt:stp-1);
        corr_sig = xcorr(erc_win_data, smarting_win_data(1:length(erc_win_data)));
        plot(corr_sig);
    end
    
    
    
   
    
    % Analyse ERC data
    oddball_analyse(EEG_erc, 1, 'ERC Data');
    
    peak_thr = 0.4;
    sel_events_smarting = extract_trigger_events(smarting_trigger_data', fs_smarting, peak_thr);
    EEG_smarting_audio_trigger = EEG_smarting;
    % The latencies of events are replaced by those extracted from audio
    % triggers
    % There are 2 triggers missing when picking up triggers from Smarting
    % audio recording. I have not looked into it closely.
    for i = 1:length(sel_events_smarting)
        if(i<length(sel_events_smarting))
            intvl(i) = sel_events_smarting(i+1) - sel_events_smarting(i);
        end
        if(intvl<1.1*fs_smarting)
            EEG_smarting_audio_trigger.event(smarting_events(i)).latency =  sel_events_smarting(i);
        end
    end
      % Analyse Smarting data
    oddball_analyse(EEG_smarting_audio_trigger, 1, 'Smarting Data Audio Triggers');
end