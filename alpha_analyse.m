function [avg_fft_open_1,avg_fft_open_2,avg_fft_close] = alpha_analyse(EEG, fs_new,dur)
    EEG = pop_resample(EEG, fs_new);
    EEG = pop_eegfiltnew(EEG, [],0.5,[],true); %HP filtering (Filter DC)
    EEG = pop_eegfiltnew(EEG, [],20); %LP filtering
    EEGfilt = EEG.data(1,:)';
    
    strt = 16;
    stp = strt+dur;
    eye_open_1_erc = EEGfilt(strt*fs_new+1:stp*fs_new);
    
    strt = stp;
    stp = strt+dur;
    eye_close_erc = EEGfilt(strt*fs_new+1:stp*fs_new);

    strt = stp;
    stp = strt+dur;
    eye_open_2_erc = EEGfilt(strt*fs_new+1:stp*fs_new);
    
    strt = 1;
    dur = 2048e-3;
    L = floor(dur*fs_new);
    stp = strt+L-1;
    
    fft_win = [];
    Hann_win = hanning(L);
    
    
    while(stp<size(eye_open_1_erc,1))
        win_data = eye_open_1_erc(strt:stp).*Hann_win;
        fft_win = [fft_win,fft(win_data)];
        strt = stp+1;
        stp =strt+L-1;
    end
    avg_fft_open_1 = mean(abs(fft_win),2);
    
    strt = 1;
    stp = strt+L-1;
    fft_win = [];
    while(stp<size(eye_open_2_erc,1))
        win_data = eye_open_2_erc(strt:stp).*Hann_win;
        fft_win = [fft_win,fft(win_data)];
        strt = stp+1;
        stp =strt+L-1;
    end
    avg_fft_open_2 = mean(abs(fft_win),2);
    
    strt = 1;
    stp = strt+L-1;
    fft_win = [];
    while(stp<size(eye_close_erc,1))
        win_data = eye_close_erc(strt:stp).*Hann_win;
        fft_win = [fft_win,fft(win_data)];
        strt = stp+1;
        stp =strt+L-1;
    end
    avg_fft_close = mean(abs(fft_win),2);        
end