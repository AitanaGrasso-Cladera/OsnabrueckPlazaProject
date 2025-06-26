# Osnabrück Plaza Project
Codes used for the Osnabrück Plaza Project developed by Aitana Grasso-Cladera, MSc. and Debora Nolte, MSc. (Neurobiopsychology Laboratory at the Institute of Cognitive Sciences; Osnabrück University), in collaboration of Aziz Muhammed Akkaya and Alina Zaidan. This project is supervised by Prof. Dr. Peter König and Prof. Dr. Tim Kietzmann.

**alignmentAndTriggerFile.m**
This script performs the aligment of EEG and eye-tracking data streams, accounting for any potential drifts in the software synchronization level.
Also, it creates the trigger file based on the results of the Face Mapper Algorithm (Pupil Lab), in order to generate a usable trigger file for EEG preprocessing and analyses.
To run this script you need the XDF file with the data streams, and different files coming from the eye-tracking data (accesible through Pupil Cloud)
The output of this script is needed in order to run OsnaPlaza_preprocessing.

**OsnaPlaza_preprocessing.m**
This script presents the preprocessing of the EEG data following the customized pipeline for the present project.
To run this script, you need the output from aligmentAndTriggerFile, the hardware trigger file, as well as the XDF file and the channel location file for the EEG system.
