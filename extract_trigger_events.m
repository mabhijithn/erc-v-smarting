function sel_events = extract_trigger_events (trigger_data, fs, thr)
    
    % Extracting energy of the signal to find peaks
    % Energy in a 60ms window extracted (the duration of a tone is 60ms)
    win_size = 60e-3;
    win_smpls = win_size*fs;
    strt = 1;
    stp  = strt+win_smpls-1;
    k = 1;
    trigger_energy = zeros(size(trigger_data,1),1);
    while(stp<size(trigger_data,1))
        win_signal = trigger_data(strt:stp);
        trigger_energy(k) = sum(win_signal.^2)/win_smpls;
        strt = strt+1;
        stp = strt+win_smpls-1;
        k = k+1;
    end



    % Using an R-peak detector to pick the peaks 
    % The function is attached in the email
    plot_r_peak_stats = 0;
  
    [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(trigger_energy,fs,plot_r_peak_stats);
    id = find(qrs_amp_raw>thr);

    % Index of oddball events in ERC data
    sel_events = qrs_i_raw(id);
end