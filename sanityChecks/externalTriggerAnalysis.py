
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import scipy.stats as stats
import mne
from pathlib import Path
import os
from collections import defaultdict
import statsmodels.api as sm

# Settings
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 500)
pd.set_option('display.precision', 15)

# Define directories (delete the path)
workingFolder = Path('')
dataFolder = Path(workingFolder, 'Data')
triggerFolder = Path(workingFolder, 'triggerFiles')
eventFolder = Path(workingFolder, 'eventFiles')

# Get participant ids
tmp = os.listdir(dataFolder) # Get the list of participants in the data folder
participants = [] # Create an empty list to store participant ids
participantFolders = [] # Create an empty list to store the paths of the participant folders
participantFiles = [] # Create an empty list to store the paths of the participant files

# Define participant ids, and paths to the folders and files
for file in tmp:
    if file.endswith('.xdf'):
        participants.append(file[:-4])
        participantFolders.append(Path(dataFolder, 'preproc_' + file[:-4]))
        participantFiles.append(Path(dataFolder, 'preproc_' + file[:-4], file[:-4] + '_with_events.set'))


# External trigger names in the EEG data
target_labels = ['A_on', 'B_on', 'C_on', 'D_on']

# A dictionary to store the data
latencies = {}

# Get the latencies from the EEG data and events files. Then calculate the residuals for external triggers.
for index, participant in enumerate(participants):
    eeg_data = mne.io.read_raw_eeglab(str(participantFiles[index]), preload=True) # Load the EEG data

    sfreq = eeg_data.info['sfreq'] # Define the sampling rate

    eeg_events, event_id = mne.events_from_annotations(eeg_data) # Get the events in the EEG data

    # Create sub-dictionaries in the main dictionary to store timestamps
    latencies[participant] = {'EEG': defaultdict(list), 'ET': defaultdict(list), 'difference': defaultdict(list)}

    # Define the drift caused by the resistor in the trigger box
    electrode_drift = 1.12195122 / 1000

    # Find the trigger events in the EEG data and create trigger files for hardware triggers
    all_values = []
    all_labels = []

    for label in target_labels:
        if label in event_id:
            event_code = event_id[label] # Find the id of the event
            matching_events = eeg_events[eeg_events[:, 2] == event_code] # Extract the relevant trigger events
            all_values.extend(matching_events[:, 0])  # add the list of numbers multiplied by the sampling frequency
            all_labels.extend(label[0] * len(matching_events[:, 0]))  # repeat the key for each number
            # Get the latencies, convert to seconds, and subtract the drift
            latencies[participant]['EEG'][label] = (matching_events[:, 0] / sfreq) - electrode_drift
        else:
            print(f"Warning: {label} not found in event_id")

    hardwareTrigger = pd.DataFrame({'latency': all_values, 'type': all_labels})
    hardwareTrigger.to_csv(Path(triggerFolder, 'hardwareTriggers_' + participant + '.csv'), index=False)


    # Import the events file
    events = pd.read_csv(Path(eventFolder, 'events_' + participant[-4:] + '.csv'))

    events = events.rename(columns={'timestamp [ns]': 'timestamps_ns'}) # Change the column name for convenience
    events['timestamps_s'] = events['timestamps_ns'] / 1e9 # Conver from nanoseconds to seconds

    events['EEG_event'] = events['name'].str.split(":").str[0] # Get the EEG event labels in the events file
    events['EEG_timestamps_s'] = events['name'].str.split(":").str[1] # Get the EEG timestamps in the events file
    events['EEG_timestamps_s'] = events['EEG_timestamps_s'].astype(float) # Change the timestamps into float

    # Define the rows to drop: first row, last to rows, and the EEG.start row
    events = events.drop([0, list(events[events['name'].str.startswith('EEG.start')].index)[0], len(events) - 2,len(events) - 1])
    events = events.reset_index(drop=True) # Drop those rows

    # Zero the starting point for both ET and EEG timestamps
    events['timestamps_s_zeroed'] = events['timestamps_s'] - events.loc[0, 'timestamps_s']
    events['EEG_timestamps_s_zeroed'] = events['EEG_timestamps_s'] - events.loc[0, 'EEG_timestamps_s']

    # Extract the trigger events in the events file
    trigger_events = events[events['EEG_event'].isin(['A', 'B', 'C', 'D'])]

    # Find the ET timestamps and calculate the residual between EEG timestamps and ET timestamps
    for label in target_labels:
        latencies[participant]['ET'][label] = trigger_events[trigger_events['EEG_event'] == label[0]]['timestamps_s_zeroed']
        latencies[participant]['difference'][label] = latencies[participant]['EEG'][label] - latencies[participant]['ET'][label]

