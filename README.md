# Yamagata_et_al_PNAS_2021_HypLinkArousalSleepHomeostasis
Codes generated for Yamagata et al, The hypothalamic link between arousal and sleep homeostasis in mice. PNAS (2021)



Please follow the title of mfiles to make figures in the paper.
(ex. For plotting the fibertip's location as shown in Fig1B, see Fig1B_Plot_fibertips.m)



Some mfiles have no figure names. These files are for below. 


%%%%%%%%    Common    %%%%%%%%

% Used in many analyses. Download and save it to your folder first. 

        resampling_new.m


% Make edf files for sleep scoring on SleepSign (you need ascii-edf converter)

        Read_EEGs_EMGs_Resample.m
        aCreateTXTfile_EEG_EMG.m (no stimulation OR stimulation timing masked)
        Read_STIM.m (readout stimulation timing)
        aCreateTXTfile_EEG_EMG_optSTIM.m (stimulation timing visualized version)
	
	
% Check SWA and EMG and check the quality of recording 	

        aPlotEEGandEMGprofilesDaysWriteOutVar.m
	

% ---------   here, Manual Scoring and FFT calculation from sleepsign  --------


% Plot EEG spectra

        ReadVSspecEEGall.m




%%%%%%%%     EMG analysis     %%%%%%%%

% Calculate EMG variance in each animals

        EMGVariance_v8_allfigs.m
	

% Calculate EMG variance for latency to arousal (Awakening delay)

        export_EMG_Variance.m
        aEMGvarianceSTIM_ty.m
	
	

% Demonstration of EMG variance distribution between vigilance states

   
   % Fig S1D
   
        extract_features_from_opto_evoked_arousals.m
	
   
   % see also 
   
        compareEMGVar_SpontaneousVsStimV8_forFun_PCAonSpectra.m




%%%%%%%%     Movie     %%%%%%%%

% Make supplimentary movie 1

        SupplMovie_videoExtraction.m

 

Δ∆
