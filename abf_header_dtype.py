import numpy as np
from abf_header_defs import *

abf_header_dtype = np.dtype([\
    ('fid_size_info', [\
        ('lFileSignature' ,  np.int32),
        ('fFileVersionNumber' ,  np.float32),
        ('nOperationMode' ,  np.int16),
        ('lActualAcqLength' ,  np.int32),
        ('nNumPointsIgnored' ,  np.int16),
        ('lActualEpisodes' ,  np.int32),
        ('lFileStartDate' ,  np.int32),
        ('lFileStartTime' ,  np.int32),
        ('lStopwatchTime' ,  np.int32),
        ('fHeaderVersionNumber' ,  np.float32),
        ('nFileType' ,  np.int16),
        ('nMSBinFormat' ,  np.int16)
        ]),
    ('f_structure', [\
        ('lDataSectionPtr' , np.int32),
        ('lTagSectionPtr' , np.int32),
        ('lNumTagEntries' , np.int32),
        ('lScopeConfigPtr' , np.int32),
        ('lNumScopes' , np.int32),
        ('_lDACFilePtr' , np.int32),
        ('_lDACFileNumEpisodes' , np.int32),
        ('sUnused001' , '|S4'),
        ('lDeltaArrayPtr' , np.int32),
        ('lNumDeltas' , np.int32),
        ('lVoiceTagPtr' , np.int32),
        ('lVoiceTagEntries' , np.int32),
        ('lUnused002' , np.int32),
        ('lSynchArrayPtr' , np.int32),
        ('lSynchArraySize' , np.int32),
        ('nDataFormat' , np.int16),
        ('nSimultaneousScan' , np.int16),
        ('lStatisticsConfigPtr' , np.int32),
        ('lAnnotationSectionPtr' , np.int32),
        ('lNumAnnotations' , np.int32),
        ('sUnused003', '|S2'),
        ('channel_count_acquired', np.int16)
        ]),
    ('trial_hierarchy', [\
        ('nADCNumChannels' , np.int16),
        ('fADCSampleInterval' , np.float32),
        ('fADCSecondSampleInterval' , np.float32),
        ('fSynchTimeUnit' , np.float32),
        ('fSecondsPerRun' , np.float32),
        ('lNumSamplesPerEpisode' , np.int32),
        ('lPreTriggerSamples' , np.int32),
        ('lEpisodesPerRun' , np.int32),
        ('lRunsPerTrial' , np.int32),
        ('lNumberOfTrials' , np.int32),
        ('nAveragingMode' , np.int16),
        ('nUndoRunCount' , np.int16),
        ('nFirstEpisodeInRun' , np.int16),
        ('fTriggerThreshold' , np.float32),
        ('nTriggerSource' , np.int16),
        ('nTriggerAction' , np.int16),
        ('nTriggerPolarity' , np.int16),
        ('fScopeOutputInterval' , np.float32),
        ('fEpisodeStartToStart' , np.float32),
        ('fRunStartToStart' , np.float32),
        ('fTrialStartToStart' , np.float32),
        ('lAverageCount' , np.int32),
        ('lClockChange' , np.int32),
        ('nAutoTriggerStrategy' , np.int16)
        ]),
    ('display_stuff',        ### 44 bytes
        [\
        ('nDrawingStrategy' , np.int16),
        ('nTiledDisplay' , np.int16),
        ('nEraseStrategy' , np.int16),
        ('nDataDisplayMode' , np.int16),
        ('lDisplayAverageUpdate' , np.int32),
        ('nChannelStatsStrategy' , np.int16),
        ('lCalculationPeriod' , np.int32),
        ('lSamplesPerTrace' , np.int32),
        ('lStartDisplayNum' , np.int32),
        ('lFinishDisplayNum' , np.int32),
        ('nMultiColor' , np.int16),
        ('nShowPNRawData' , np.int16),
        ('fStatisticsPeriod' , np.float32),
        (' lStatisticsMeasurements' , np.int32),
        ('nStatisticsSaveStrategy' , np.int16)
        ]),
    ('hardware_inf',            ### 16 bytes
        [\
        ('fADCRange' , np.float32),
        ('fDACRange' , np.float32),
        ('lADCResolution' , np.int32),
        ('lDACResolution' , np.int32)
        ]),
    ('environment_inf',         ### 118 bytes
        [\
        ('nExperimentType' , np.int16),
        ('_nAutosampleEnable' , np.int16),
        ('_nAutosampleADCNum' , np.int16),
        ('_nAutosampleInstrument' , np.int16),
        ('_fAutosampleAdditGain' , np.float32),
        ('_fAutosampleFilter' , np.float32),
        ('_fAutosampleMembraneCap' , np.float32),
        ('nManualInfoStrategy' , np.int16),
        ('fCellID1' , np.float32),
        ('fCellID2' , np.float32),
        ('fCellID3' , np.float32),
        ('sCreatorInfo',  np.str_, ABF_CREATORINFOLEN), # is cool for strings
        ('_sFileComment', np.str_, ABF_OLDFILECOMMENTLEN),
        ('nFileStartMillisecs' , np.int16),
        ('nCommentsEnable' , np.int16),
        ('sUnused003a', np.str_, 8)
        ]),
    ('multi-chan_inf',          ###  (1044 (160 + 384 + 488 + 12) bytes)
        [\
        ('nADCPtoLChannelMap', np.int16, ABF_ADCCOUNT),
        ('nADCSamplingSeq', np.int16, ABF_ADCCOUNT),
        ('sADCChannelName', ('|S' + str(ABF_ADCNAMELEN), ABF_ADCCOUNT) ),  ####THIS IS the way to do arrays of strings
        ('sADCUnits', ('|S' + str(ABF_ADCUNITLEN) , ABF_ADCCOUNT)),
        ('fADCProgrammableGain', np.float32, ABF_ADCCOUNT),
        ('fADCDisplayAmplification', np.float32, ABF_ADCCOUNT),
        ('fADCDisplayOffset'      , np.float32, ABF_ADCCOUNT),
        ('fInstrumentScaleFactor' , np.float32, ABF_ADCCOUNT),
        ('fInstrumentOffset'      , np.float32, ABF_ADCCOUNT),
        ('fSignalGain', np.float32, ABF_ADCCOUNT),
        ('fSignalOffset', np.float32, ABF_ADCCOUNT),
        ('fSignalLowpassFilter', np.float32, ABF_ADCCOUNT),
        ('fSignalHighpassFilter', np.float32, ABF_ADCCOUNT),
        ('sDACChannelName', ('|S' + str(ABF_DACNAMELEN),ABF_DACCOUNT)),
        ('sDACChannelUnits', ('|S' + str(ABF_DACUNITLEN),ABF_DACCOUNT)),
        ('fDACScaleFactor', np.float32, ABF_DACCOUNT),
        ('fDACHoldingLevel', np.float32, ABF_DACCOUNT),
        ('nSignalType', np.int16),
        ('sUnused004', np.str_, 10)
        ]),
    ('synch_timer_outs',         ### (14 bytes)
        [\
        ('nOUTEnable' , np.int16),
        ('nSampleNumberOUT1' , np.int16),
        ('nSampleNumberOUT2' , np.int16),
        ('nFirstEpisodeOUT' , np.int16),
        ('nLastEpisodeOUT' , np.int16),
        ('nPulseSamplesOUT1' , np.int16),
        ('nPulseSamplesOUT2' , np.int16)
        ]),
    ('epoch_waveform_pulses',    ### (184 bytes)
        [\
        ('nDigitalEnable' , np.int16),
        ('_nWaveformSource' , np.int16),
        ('nActiveDACChannel' , np.int16),
        ('_nInterEpisodeLevel' , np.int16),
        ('_nEpochType', np.int16, ABF_EPOCHCOUNT),
        ('_fEpochInitLevel', np.float32, ABF_EPOCHCOUNT),
        ('_fEpochLevelInc', np.float32, ABF_EPOCHCOUNT),
        ('_nEpochInitDuration', np.int16, ABF_EPOCHCOUNT),
        ('_nEpochDurationInc', np.int16, ABF_EPOCHCOUNT),
        ('nDigitalHolding' , np.int16),
        ('nDigitalInterEpisode' , np.int16),
        ('nDigitalValue', np.int16, ABF_EPOCHCOUNT),
        ('sUnavailable1608', np.str_, 4),
        ('nDigitalDACChannel' , np.int16),
        ('sUnused005', np.str_, 6)
        ]),
    ('DAC_output_file',          ### (98 bytes)
        [\
        ('_fDACFileScale' , np.float32),
        ('_fDACFileOffset' , np.float32),
        ('sUnused006', np.str_, 2),
        ('_nDACFileEpisodeNum' , np.int16),
        ('_nDACFileADCNum' , np.int16),
        ('_sDACFilePath', np.str_, ABF_DACFILEPATHLEN)
        ]),
    ('presweep_pulse_train',     ### (44 bytes)
        [\
        ('_nConditEnable' , np.int16),
        ('_nConditChannel' , np.int16),
        (' _lConditNumPulses' , np.int32),
        ('_fBaselineDuration' , np.float32),
        ('_fBaselineLevel' , np.float32),
        ('_fStepDuration' , np.float32),
        ('_fStepLevel' , np.float32),
        ('_fPostTrainPeriod' , np.float32),
        ('_fPostTrainLevel' , np.float32),
        ('sUnused007', np.str_, 12)
        ]),
    ('variable_param_user_list', ### ( 82 bytes)
        [\
        ('_nParamToVary' , np.int16),
        ('_sParamValueList', np.str_, ABF_VARPARAMLISTLEN)
        ]),
    ('autopeak_msr',             ### (36 bytes)
        [\
        ('_nAutopeakEnable' , np.int16),
        ('_nAutopeakPolarity' , np.int16),
        ('_nAutopeakADCNum' , np.int16),
        ('_nAutopeakSearchMode' , np.int16),
        (' _lAutopeakStart' , np.int32),
        (' _lAutopeakEnd' , np.int32),
        ('_nAutopeakSmoothing' , np.int16),
        ('_nAutopeakBaseline' , np.int16),
        ('_nAutopeakAverage' , np.int16),
        ('sUnavailable1866', np.str_, 2),
        (' _lAutopeakBaselineStart' , np.int32),
        (' _lAutopeakBaselineEnd' , np.int32),
        (' _lAutopeakMeasurements' , np.int32)
        ]),
    ('channel_math',             ### (52 bytes)
        [\
        ('nArithmeticEnable' , np.int16),
        ('fArithmeticUpperLimit' , np.float32),
        ('fArithmeticLowerLimit' , np.float32),
        ('nArithmeticADCNumA' , np.int16),
        ('nArithmeticADCNumB' , np.int16),
        ('fArithmeticK1' , np.float32),
        ('fArithmeticK2' , np.float32),
        ('fArithmeticK3' , np.float32),
        ('fArithmeticK4' , np.float32),
        ('sArithmeticOperator', np.str_, ABF_ARITHMETICOPLEN),
        ('sArithmeticUnits', np.str_, ABF_ARITHMETICUNITSLEN),
        ('fArithmeticK5' , np.float32),
        ('fArithmeticK6' , np.float32),
        ('nArithmeticExpression' , np.int16),
        ('sUnused008', np.str_, 2)
        ]),
     ('on-line_sub',     ### (34 bytes)
        [\
        ('_nPNEnable' , np.int16),
        ('nPNPosition' , np.int16),
        ('_nPNPolarity' , np.int16),
        ('nPNNumPulses' , np.int16),
        ('_nPNADCNum' , np.int16),
        ('_fPNHoldingLevel' , np.float32),
        ('fPNSettlingTime' , np.float32),
        ('fPNInterpulse' , np.float32),
        ('sUnused009', np.str_, 12)
        ]),
     ('unused_space',            ### (54 bytes)
        [\
        ('_nListEnable' , np.int16),
        ('nBellEnable', np.int16, ABF_BELLCOUNT),
        ('nBellLocation', np.int16, ABF_BELLCOUNT),
        ('nBellRepetitions', np.int16, ABF_BELLCOUNT),
        ('nLevelHysteresis' , np.int16),
        ('lTimeHysteresis' , np.int32),
        ('nAllowExternalTags' , np.int16),
        ('nLowpassFilterType', np.str_, ABF_ADCCOUNT),
        ('nHighpassFilterType', np.str_, ABF_ADCCOUNT),
        ('nAverageAlgorithm' , np.int16),
        ('fAverageWeighting' , np.float32),
        ('nUndoPromptStrategy' , np.int16),
        ('nTrialTriggerSource' , np.int16),
        ('nStatisticsDisplayStrategy' , np.int16),
        ('nExternalTagType' , np.int16),
        ('lHeaderSize' , np.int32),
        ('FileDuration' , np.float64), ######################not sure about type
        ('nStatisticsClearStrategy' , np.int16)
        ]),
     ('ext_f_structure',         ### (8 * 2 = 16 + 100 bytes)
        [\
        ('lDACFilePtr', np.int32, ABF_WAVEFORMCOUNT),
        ('lDACFileNumEpisodes', np.int32, ABF_WAVEFORMCOUNT),
        ('sUnused010', np.str_, 10)
        ]),
     ('ext_multi-chan_inf',      ### (8 * 4 = 32 + 100 bytes)
        [\
        ('fDACCalibrationFactor', np.float32, ABF_DACCOUNT),
        ('fDACCalibrationOffset', np.float32, ABF_DACCOUNT),
        ('sUnused011', np.str_, 30),
        ('lEpochPulsePeriod', np.int32, (ABF_WAVEFORMCOUNT, ABF_EPOCHCOUNT)), ## 4 * 2 * 10 = 80 bytes
        ('lEpochPulseWidth' , np.int32, (ABF_WAVEFORMCOUNT, ABF_EPOCHCOUNT)), ## 4 * 2 * 10 = 80 bytes
        ]),
     ('ext_epoch_waveform_pulses',   ### ( 186 * 2 = 368 + 100 bytes)
        [\
        ('nWaveformEnable', np.int16, ABF_WAVEFORMCOUNT),
        ('nWaveformSource', np.int16, ABF_WAVEFORMCOUNT),
        ('nInterEpisodeLevel', np.int16, ABF_WAVEFORMCOUNT),
        ('nEpochType', np.int16, (ABF_WAVEFORMCOUNT, ABF_EPOCHCOUNT)),
        ('fEpochInitLevel', np.float32, (ABF_WAVEFORMCOUNT, ABF_EPOCHCOUNT)),
        ('fEpochLevelInc', np.float32, (ABF_WAVEFORMCOUNT, ABF_EPOCHCOUNT)),
        ('lEpochInitDuration', np.int32, (ABF_WAVEFORMCOUNT, ABF_EPOCHCOUNT)),
        ('lEpochDurationInc', np.int32, (ABF_WAVEFORMCOUNT, ABF_EPOCHCOUNT)),
        ('nDigitalTrainValue', np.int16, ABF_EPOCHCOUNT), ## 2 * 10 = 20 bytes
        ('nDigitalTrainActiveLogic', np.int16),   ## 2 bytes
        ('sUnused012', np.str_, 18)
        ]),
    ('ext_DAC_output_file',     ### (270 * 2 = 540 + 12 bytes)
        [\
        ('fDACFileScale', np.float32, ABF_WAVEFORMCOUNT),
        ('fDACFileOffset', np.float32, ABF_WAVEFORMCOUNT),
        ('lDACFileEpisodeNum', np.int32, ABF_WAVEFORMCOUNT),
        ('nDACFileADCNum', np.int16, ABF_WAVEFORMCOUNT),
        ('sDACFilePath', ('|S' + str(ABF_PATHLEN),ABF_WAVEFORMCOUNT)),
        ('sUnused013', np.str_, 12)
        ]),
    ('ext_presweep_pulse_train', ### (2 * 32 = 64 + 36 bytes)
        [\
        ('nConditEnable', np.int16, ABF_WAVEFORMCOUNT),
        ('lConditNumPulses', np.int32, ABF_WAVEFORMCOUNT),
        ('fBaselineDuration', np.float32, ABF_WAVEFORMCOUNT),
        ('fBaselineLevel', np.float32, ABF_WAVEFORMCOUNT),
        ('fStepDuration', np.float32, ABF_WAVEFORMCOUNT),
        ('fStepLevel', np.float32, ABF_WAVEFORMCOUNT),
        ('fPostTrainPeriod', np.float32, ABF_WAVEFORMCOUNT),
        ('fPostTrainLevel', np.float32, ABF_WAVEFORMCOUNT),
        ('sUnused014', np.str_, 40)
        ]),
    ('ext_variable_param_user_list', ### (260 * 4 = 1040 + 56 bytes)
        [\
        ('nULEnable', np.int16, ABF_USERLISTCOUNT),
        ('nULParamToVary', np.int16, ABF_USERLISTCOUNT),
        ('sULParamValueList', ('|S' + str(ABF_USERLISTLEN),ABF_USERLISTCOUNT)),
        ('nULRepeat', np.int16, ABF_USERLISTCOUNT),
        ('sUnused015', np.str_, 48)
        ]),
    ('ext_on-line_sub',  ### (10 * 2 = 20 + 36 bytes)
        [\
        ('nPNEnable', np.int16, ABF_WAVEFORMCOUNT),
        ('nPNPolarity', np.int16, ABF_WAVEFORMCOUNT),
        ('nPNADCNum', np.int16, ABF_WAVEFORMCOUNT),
        ('fPNHoldingLevel', np.float32, ABF_WAVEFORMCOUNT),
        ('sUnused016', np.str_, 36)
        ]),
    ('ext_environment_inf',  ### ( (16 * 20) + 128 = 448 + 128 bytes)
        [\
        ('nTelegraphEnable', np.int16, ABF_ADCCOUNT),
        ('nTelegraphInstrument', np.int16, ABF_ADCCOUNT),
        ('fTelegraphAdditGain', np.float32, ABF_ADCCOUNT),
        ('fTelegraphFilter', np.float32, ABF_ADCCOUNT),
        ('fTelegraphMembraneCap', np.float32, ABF_ADCCOUNT),
        ('nTelegraphMode', np.int16, ABF_ADCCOUNT),
        ('nManualTelegraphStrategy', np.int16, ABF_ADCCOUNT),
        ('nAutoAnalyseEnable', np.int16),
        ('sAutoAnalysisMacroName', np.str_, ABF_MACRONAMELEN),
        ('sProtocolPath', np.str_, ABF_PATHLEN),
        ('sFileComment', np.str_, ABF_FILECOMMENTLEN),
        ('sUnused017', np.str_, 128)
        ]),
    ('ext_stats_msr', ### ( 216 + 204 = 420 bytes)
        [\
        ('nStatsEnable' , np.int16),
        ('nStatsActiveChannels', np.uint16), ############not sure about type
        ('nStatsSearchRegionFlags', np.uint16),
        ('nStatsSelectedRegion' , np.int16),
        ('_nStatsSearchMode' , np.int16),
        ('nStatsSmoothing' , np.int16),
        ('nStatsSmoothingEnable' , np.int16),
        ('nStatsBaseline' , np.int16),
        ('lStatsBaselineStart' , np.int32),
        ('lStatsBaselineEnd' , np.int32),
        ('lStatsMeasurements', np.int32, ABF_STATS_REGIONS),
        ('lStatsStart', np.int32, ABF_STATS_REGIONS),
        ('lStatsEnd', np.int32, ABF_STATS_REGIONS),
        ('nRiseBottomPercentile', np.int16, ABF_STATS_REGIONS),
        ('nRiseTopPercentile', np.int16, ABF_STATS_REGIONS),
        ('nDecayBottomPercentile', np.int16, ABF_STATS_REGIONS),
        ('nDecayTopPercentile', np.int16, ABF_STATS_REGIONS),
        ('nStatsChannelPolarity', np.int16, ABF_ADCCOUNT),
        ('nStatsSearchMode', np.int16, ABF_STATS_REGIONS),
        ('sUnused018', np.str_, 156)
        ]),
    ('ext_app_ver_data',   ### ( 2 * 4 + 8 = 16 bytes )
        [\
        ('nMajorVersion' , np.int16),
        ('nMinorVersion' , np.int16),
        ('nBugfixVersion' , np.int16),
        ('nBuildVersion' , np.int16),
        ('sUnused019', np.str_, 8)
        ]),
    ('ext_LTP_prot', ### ( 2 + ( 4 * 2 ) + 4 = 14 bytes )
        [\
        ('nLTPType' , np.int16),
        ('nLTPUsageOfDAC', np.int16, ABF_WAVEFORMCOUNT),
        ('nLTPPresynapticPulses', np.int16, ABF_WAVEFORMCOUNT),
        ('sUnused020', np.str_, 4)
        ]),
    ('ext_digiD_trig_out_flag', ### ( 2 + 6 = 8 bytes )
        [\
        ('nDD132xTriggerOut' , np.int16),
        ('sUnused021', np.str_, 6)
        ]),
    ('ext_epoch_resis', ### ( ( 2 * 10 ) + ( 2 * 2 ) + 32 = 56 bytes )
        [\
        ('sEpochResistanceSignalName', ('|S' + str(ABF_ADCNAMELEN),ABF_WAVEFORMCOUNT)),
        ('nEpochResistanceState', np.int16, ABF_WAVEFORMCOUNT),
        ('sUnused022', np.str_, 16)
        ]),
    ('ext_alt_episod_mode', ### ( ( 2 * 10 * 2 ) + 2 + 32 = 104 bytes )
        [\
        ('nAlternateDACOutputState' , np.int16),
        ('nAlternateDigitalValue', np.int16, ABF_EPOCHCOUNT),
        ('nAlternateDigitalTrainValue', np.int16, ABF_EPOCHCOUNT),
        ('nAlternateDigitalOutputState' , np.int16),
        ('sUnused023', np.str_, 14)
        ]),
    ('ext_post-procss', ### ( (4*16) + 16 = 80)
        [\
        ('fPostProcessLowpassFilter', np.float32, ABF_ADCCOUNT),
        ('nPostProcessLowpassFilterType', np.str_, ABF_ADCCOUNT)
        ]),
    ('ext_unused', [('sUnused2048', np.str_, 130)]) ####### is 130
])
