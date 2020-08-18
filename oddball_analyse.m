function oddball_analyse(EEG, scaling, title_str)
    
% Select only relevent channels of data
    EEG = pop_select( EEG,'channel',{'Channel 1'});
    EEG = pop_eegfiltnew(EEG, [],0.5,[],true); %HP filtering
    EEG = pop_eegfiltnew(EEG, [],20); %LP filtering

    EEGB = pop_epoch( EEG, {  '33025'  }, [-0.2  1.0], 'newname', 'GDF file epochs', 'epochinfo', 'yes'); %025 Extract Background
    EEGD = pop_epoch( EEG, {  '33026'  }, [-0.2  1.0], 'newname', 'GDF file epochs', 'epochinfo', 'yes'); %026 Extract Target High tone

    EEGB = pop_rmbase( EEGB, [-200   0]);
    EEGD = pop_rmbase( EEGD, [-200   0]);

    x_axis = EEG.srate.*(0.2+[-0.2:0.1:0.5,0.8,1.0]);

    figure
    for pl = 1:1
        subplot(1,1,pl)
        plot(scaling.*mean(EEGD.data(pl,:,:),3),'r'); hold on;
        plot(scaling.*mean(EEGB.data(pl,:,:),3),'b'); %hold on;
        title(sprintf('Channel %d of %s', pl, title_str))
        xlabel('ms');
        ylabel('Average across epochs');
        legend('Target (high)','Background')
        set(gca, 'XTick', x_axis);
        set(gca,'XTickLabel',{'-200'; '-100';'S';'100';'200';'300';'400';'500';'800';'1000'});
    end
    
end