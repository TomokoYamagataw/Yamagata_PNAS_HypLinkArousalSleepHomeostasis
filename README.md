# Yamagata_et_al_PNAS_2021_HypLinkArousalSleepHomeostasis
Codes generated for Yamagata et al, The hypothalamic link between arousal and sleep homeostasis in mice. PNAS (2021)



Please follow the title of mfiles to make figures in the paper.
(ex. For plotting the fibertip's location as shown in Fig1B, see Fig1B_Plot_fibertips.m)



Some basic mfiles have no figure names. These files are for below. 


% Used in many analyses. Download and save it to your folder first. 
        resampling_new.m

% Make edf files for scoring on SleepSign (additionally, you need ascii to edf file converter)
        Read_EEGs_EMGs_Resample.m
        aCreateTXTfile_EEG_EMG.m (no stimulation OR stimulation timing masked)
        Read_STIM.m (readout stimulation timing)
        aCreateTXTfile_EEG_EMG_optSTIM.m (stimulation timing visualized version)
	
% Check SWA and EMG before scoring	
        aPlotEEGandEMGprofilesDaysWriteOutVar.m

% Calculate EMG variance for wake latency analysis
        export_EMG_Variance.m
        aEMGvarianceSTIM_ty.m
	
% Make supplimentary movie 1
        SupplMovie_videoExtraction.m


% Calculate EMG variance in each animals
        EMGVariance_v8_allfigs.m

% Demonstration of EMG variance distribution between vigilance states
   % Fig S1B
        extract_features_from_opto_evoked_arousals.m
   % see also 
        compareEMGVar_SpontaneousVsStimV8_forFun_PCAonSpectra.m
  
% Latency to Arousal (Awakening delay)
   % Fig 1H & Fig S1C, LPO vs nonLPO (NREM vs REM)
        LPOvsnonLPO_fig1.m
	
   % Fig 2D: HSP vs LSP
        highLowSleepPressureTArousal_fig2.m
	
   % Fig S3F: Dex sedation
        sedationTArousal_fig2.m
	
   % Fig 5A: 1,2,5,10,20Hz,sham
        compare1nevs5vs10Hz_V5.m
 

Δ∆
