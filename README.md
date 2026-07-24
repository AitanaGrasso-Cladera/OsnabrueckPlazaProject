# Osnabrück Plaza Project
Codes used for the Osnabrück Plaza Project developed by Aitana Grasso-Cladera, MSc. and Debora Nolte, MSc. (Neurobiopsychology Laboratory at the Institute of Cognitive Sciences; Osnabrück University), in collaboration of Aziz Muhammed Akkaya and Alina Zaidan. This project is supervised by Prof. Dr. Peter König and Prof. Dr. Tim Kietzmann. 

The data used in the project can be found under: https://osf.io/dvgue/overview

**alignmentAndTriggerFile.m**
This script performs the aligment of EEG and eye-tracking data streams, accounting for any potential drifts in the software synchronization level.
Also, it creates the trigger file based on the results of the Face Mapper Algorithm (Pupil Lab), in order to generate a usable trigger file for EEG preprocessing and analyses.
To run this script you need the XDF file with the data streams, and different files coming from the eye-tracking data (accesible through Pupil Cloud).
The output of this script is needed in order to run OsnaPlaza_preprocessing and position event markers.

**findingHardwareTriggers.m**
This script allows the identification of the artifact triggers generated in the EEG signal, by using temporal information from the artifacts on the eye-tracking video.
After identifying the artifacts, it generates a trigger file with the events corresponding to the different rounds of hardware trigger during the experimental session.
To run this script you need the XDF file with the data streams, and the event file coming from the eye-tracking data (accesible through Pupil Cloud).
The output of this script is needed in order to run OsnaPlaza_preprocessing and clean the data.

**OsnaPlaza_preprocessing.m**
This script presents the preprocessing of the EEG data following the customized pipeline for the present project.
To run this script, you need the output from aligmentAndTriggerFile, the hardware trigger file, as well as the XDF file and the channel location file for the EEG system.

The rest of the scripts should be used in the following order:
* 1. _Saccade and Blink Dynamics During Free Exploration_ section and plot _Figure 3_
- a. figure3_SaccBlinkDurations.m
2. _Contrasting Saccade- and Blink-Offset ERPs_ section and plot _Figure 4_
  a. epochNoUnfold.m
  b. ERPimageSortedBlinkDur.m
  c. ERPimageSorderSaccDur.m
  d. ERPs_Unfold_OSPlaza.m
  e. ERPs_Unfold_OSPlaza_TF.m
  f. figure4_SaccBlinkSorted.m
3. _Comparing Saccade- and Blink-Locked ERPs_ section and plot _Figure 5_
  a. averageSessions.m
  b. rootToMeanRatio.m.
  c. figure5_6_A_Activity.m (requires the function compute_CI.m)
  d. figure5_6_BandC_CosineSimilarityMatrix.m
  e. figure5_6_D_TheoreticalMatrix.m
  f. theoreticalTemplateCorrelation.m (this script also generates Panel E of Figures 5 and 6).
  g. sign_flippedCorrelation.m
