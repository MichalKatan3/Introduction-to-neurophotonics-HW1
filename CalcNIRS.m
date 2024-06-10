function [dHbR , dHbO, fig] = CalcNIRS(dataFile, SDS, tissueType, plotChannelIdx, extinctionCoefficientsFile, DPFperTissueFile, relDPFfile )
for ind = 1:length(plotChannelIdx)
    Lambda1 = dataFile.SD.Lambda(1); % [nm]
    Lambda2 = dataFile.SD.Lambda(2); % [nm]
    Intensity_1 = dataFile.d(:,1:20);  % 20 channels for each wavelength
    Intensity_2 = dataFile.d(:,21:40); 

    DPF = DPFperTissueFile.DPF(tissueType); % [#]
    relative_DPF_index_1 = find(relDPFfile.wavelength == Lambda1);
    relative_DPF_index_2 = find(relDPFfile.wavelength == Lambda2);
    relative_DPF1 = relDPFfile.relDPFcoeff(relative_DPF_index_1);
    relative_DPF2 = relDPFfile.relDPFcoeff(relative_DPF_index_2);

    ch = plotChannelIdx(ind); % the current channel of plotChannelIdx
    OD1 = log10(Intensity_1(1,ch)./Intensity_1(:,ch));
    OD2 = log10(Intensity_2(1,ch)./Intensity_2(:,ch));
    OD = [OD1'; OD2'];
    Leff1 = DPF*relative_DPF1*SDS; % cm
    Leff2 = DPF*relative_DPF2*SDS; % cm

    row_idx1 = find(extinctionCoefficientsFile.wavelength == Lambda1);
    row_idx2 = find(extinctionCoefficientsFile.wavelength == Lambda2);

    ext_HbO_1 = extinctionCoefficientsFile.HbO2(row_idx1);
    ext_HbR_1 = extinctionCoefficientsFile.HHb(row_idx1);
    ext_HbO_2 = extinctionCoefficientsFile.HbO2(row_idx2);
    ext_HbR_2 = extinctionCoefficientsFile.HHb(row_idx2);

    ext = [ext_HbR_1, ext_HbO_1; ext_HbR_2, ext_HbO_2]/log(10); 

    del_conc = ext^(-1)* OD / Leff1; 

    dHbR = del_conc(1,:);
    dHbO = del_conc(2,:);

    %% FFT
    data = Intensity_1(:,ch);
    N = length(data);
    dt = dataFile.t(2)-dataFile.t(1);
    df = (1/N)*(1./dt);    
    f = 0:df:(df*length(data)-df);  

    %fft:
    fft_intensity = abs(fft(data)/N);

    %Heart rate (HR):
    window_HR = zeros(1,length(f));
    window_HR(f>1 & f<2)=1; % The "normal" range of HR (60-120 bpm equals to 1-2 Hz)
    [HR_freq_val, HR_freq_idx] = max(window_HR' .* fft_intensity); % find the HR peak

    % SNR:
    f_max = length(f)/2;
    noise = mean(fft_intensity(f > 2.5 & f<f_max));
    signal = HR_freq_val*2;

    SNR{ind,1} = round(signal/noise,2);


    %% Plottings
    darkRed_color = [120/255, 0, 0];
    blue_color = [0/255, 48/255, 150/255];
    figWidth = 25;
    figHeight = 15;
    tissue_type = DPFperTissueFile.Tissue{tissueType,1}; % for the title
    tissue_type(tissue_type=='_')=' ';
    
    fig = figure;
    set(fig, 'Units', 'centimeters');    
    set(fig, 'Position', [5, 5, figWidth, figHeight]);
    sgtitle({['Channel ', num2str(plotChannelIdx(ind))]; ['Tissue type: ', tissue_type]}, 'FontWeight', 'bold')
    subplot(2,1,1)
        plot(dataFile.t,del_conc(1,:), 'color', blue_color)
        hold on
        plot(dataFile.t,del_conc(2,:), 'color', darkRed_color)
        legend('HbR', 'HbO', 'Location', 'northeastoutside')
        xlabel('Time [sec]')
        ylabel('\Delta Hb')
        title('\Delta Hb over time')
        xlim([0,dataFile.t(end)])
    subplot(2,1,2)
        semilogy(f(2:f_max), fft_intensity(2:f_max), 'color', 'k')
        hold on
        plot(df*HR_freq_idx, HR_freq_val, 'ro', 'LineWidth', 2)
        xlabel('Frequency [Hz]')
        title('Fourier domain')        
        xlim([0,df*f_max])
        text(f_max*df*0.72, max(fft_intensity(2:f_max))*0.1, ['SNR = ', num2str(SNR{ind,1})], 'color', darkRed_color, 'FontWeight', 'bold', 'FontSize', 14);

end
end