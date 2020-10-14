function plot_fft_analysis(fftstoplot, savefig, savepath)
    if(nargin>1)
        if(savefig)
            if(isempty(savepath))
                mkdir('tmpfig');
                savepath = 'tmpfig';
            end
        end
    else
        savefig = 0;
    end
    
    width = 3.25;
    height=2.75;
    
    figure('Units','inches',...
    'Position',[0 0 width height],...
    'PaperPositionMode','auto');
    for i = 1:size(fftstoplot,2)        
        L = fftstoplot(i).L;
        L_by_2 = floor(L/2);
        fftdata = fftstoplot(i).data;
        title_str = fftstoplot(i).title;
        fs = fftstoplot(i).fs;
        savestr = fftstoplot(i).savestr;

        %calculate frequency bins with FFT
        df=fs/L; %frequency resolution
        sampleIndex = 0:L-1; %raw index for FFT plot
        f=sampleIndex*df; %x-axis index converted to frequencies
        plot(f(1:L_by_2),abs(fftdata(1:L_by_2)));
        ylabel('FFT Magnitude (Avg)');
        xlabel('Frequency (Hz)');
        axis([0 40 0 170]);
        set(gca,'XTick',[0 5 10 15 20 40]);
        grid on;
        hold on;                
    end
    legend({fftstoplot(:).title});
    if(savefig)
        saveas(gcf, fullfile(savepath,savestr),'epsc');
    end
end