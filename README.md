# Yamagata_et_al_PNAS_2021_HypLinkArousalSleepHomeostasis
Codes generated for Yamagata et al, The hypothalamic link between arousal and sleep homeostasis in mice. PNAS (2021)



Please follow the title of mfiles to make figures in the paper.
(ex. For plotting the fibertip's location as shown in Fig1B, see Fig1B_Plot_fibertips.m)



Some basic mfiles have no figure names. These files are for below. 


% Used in many analyses. Download and save it to your folder first. 

        resampling_new.m


% Make edf files for sleep scoring on SleepSign (you need ascii-edf converter)

        Read_EEGs_EMGs_Resample.m
        aCreateTXTfile_EEG_EMG.m (no stimulation OR stimulation timing masked)
        Read_STIM.m (readout stimulation timing)
        aCreateTXTfile_EEG_EMG_optSTIM.m (stimulation timing visualized version)
	
	
% Check SWA and EMG and check the quality of recording 	

        aPlotEEGandEMGprofilesDaysWriteOutVar.m
	

%%%%%%     Manual Scoring and FFT calculation from sleepsign   %%%%%%%%


% Calculate EEG spectra

        ReadVSspecEEGall.m


% Calculate EMG variance for wake latency analysis

        export_EMG_Variance.m
        aEMGvarianceSTIM_ty.m
	
	
% Make supplimentary movie 1

        SupplMovie_videoExtraction.m


% Calculate EMG variance in each animals

        EMGVariance_v8_allfigs.m
	

% Demonstration of EMG variance distribution between vigilance states

   
   % Fig S1D
   
        extract_features_from_opto_evoked_arousals.m
	
   
   % see also 
   
        compareEMGVar_SpontaneousVsStimV8_forFun_PCAonSpectra.m


  
% Latency to Arousal (Awakening delay)

   
   % Fig 1H, J, Fig S1E, LPO vs nonLPO (NREM vs REM)
   
        LPOvsnonLPO_fig1.m
	
	
  
   % Fig 3G: HSP vs LSP
   
        highLowSleepPressureTArousal_fig2.m
	
	
   
   % Fig S2: 1,2,5,10,20Hz,sham
   
        compare1nevs5vs10Hz_V5.m
	
   
   % Fig S3G: Dex sedation
   
        sedationTArousal_fig2.m
	

 
 

Δ∆
