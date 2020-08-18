clc;
clear;
close all;
if 1
    % 
    % Location on your PC of the data (Feb-05)
    main_fldr = '/home/abhijith/MEGAsync/PhD/2019/masters_thesis_micas/ERC_vs_Smarting/Feb-12';
    Fs = 120; %1000;
    fl = 20;
    fh = 0.5;
    [ BP_equirip ] = cnstr_bpfilter(Fs, fl, fh);

    % Add the peak-detection function to the MATLAB path
    % The toolbox is attached with the email
    addpath('/home/abhijith/Documents/MATLAB/toolboxes/pan_tompkin');
    
    % Following files are with phones on airplance mode
    oddball_erc_file = 'oddballTest1_2020_02_15.csv';
    oddball_smarting_file = 'ercboard-v-smarting-oddball-test-1-[2020.02.12-15.21.46].gdf';
   
   
    EEG_smarting = pop_biosig(fullfile(main_fldr, oddball_smarting_file));
    EEG_ercboard = csvread(fullfile(main_fldr, oddball_erc_file),2,0);
    EEG_erc_data = EEG_ercboard(:,2);
    
    fs_smarting = EEG_smarting.srate;
    fs_erc = 5000;
    fs_new = 500;
    
    % Prepare EEGLAB structure from ERC data
    EEG_erc = EEG_smarting;
    EEG_erc.data = EEG_ercboard(:,2)';
    EEG_erc.srate = fs_erc;
    EEG_erc.event = EEG_smarting.event;
    EEG_erc.pnts = size(EEG_erc.data,1);
    EEG_erc.times = (0:EEG_erc.pnts-1).*(1/EEG_erc.srate);
    EEG = pop_select( EEG_erc,'channel',{'Channel 1'});
    EEG = pop_eegfiltnew(EEG_erc, [],0.5,[],true); %HP filtering
    EEG = pop_eegfiltnew(EEG_erc, [],20); %LP filtering
    
%     EEG_erc_data_filtered = filtfilt(BP_equirip.numerator, 1, double(EEG_erc_data));
    

    
    smarting_data = EEG_smarting.data(1,:); % EEG data
    smarting_data_filtered = filtfilt(BP_equirip.numerator, 1, double(smarting_data'));
    
    dur = 2048e-3;
    L = floor(dur*fs_new);
    Hann_window = hann(L);
    
    strt = 1;
    stp = strt+L-1;
    %calculate frequency bins with FFT
    df=fs_new/L; %frequency resolution
    sampleIndex = 0:L-1; %raw index for FFT plot
    f=sampleIndex*df; %x-axis index converted to frequencies
    while stp<size(smarting_data_filtered,1)
        win_data = smarting_data_filtered(strt_stp)*Hann_window;
        win_fft = fft(win_data);
        alpha_power = win_fft(f>8 && f<12);
    end
end
    