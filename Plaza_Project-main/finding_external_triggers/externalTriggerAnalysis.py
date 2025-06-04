
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

# Define directories
workingFolder = Path('/Users/azizakkaya/Library/CloudStorage/GoogleDrive-aziz.muhammed64@gmail.com/My Drive/ERASMUS-Osnabruck/lab_files/plaza_project/Plaza_Project-main/finding_external_triggers')
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



### TRIGGER EVENT PARAMETERS ###
### Create an empty data frame to store external trigger parameters like mean and median
df = pd.DataFrame()

# Define the column names
df[['participant', 'target label', 'mean', 'median']] = None

i = 0 # Get a counter

# Calculate the parameters
for index, participant in enumerate(participants):
    for label in target_labels:
        df.loc[i, 'participant'] = participant
        df.loc[i, 'target label'] = label
        df.loc[i, 'mean'] = np.mean(latencies[participant]['difference'][label])
        df.loc[i, 'median'] = np.median(latencies[participant]['difference'][label])
        i = i + 1

# Store the data frame as a csv file
df.to_csv(Path(workingFolder, 'hardware_synch_parameters_with_electrode_drift.csv'))

### PLOTS ###

# Create a set of plots to show residuals
fig, axs = plt.subplots(2, 2, figsize=(10, 8)) # Define the plot plane

lines = {}

for participant in participants:
    for label in target_labels:
        x = list(range(len(latencies[participant]['difference'][label])))
        y = latencies[participant]['difference'][label]
        if label == 'A_on':
            line = axs[0, 0].plot(x, y, label=participant)
        elif label == 'B_on':
            line = axs[0, 1].plot(x, y, label=participant)
        elif label == 'C_on':
            line = axs[1, 0].plot(x, y, label=participant)
        else:
            line = axs[1, 1].plot(x, y, label=participant)

        # Store the first line for each participant (they should all be the same color)
        if participant not in lines:
            lines[participant] = line[0]


axs[0, 0].set_title('A')
axs[0, 1].set_title('B')
axs[1, 0].set_title('C')
axs[1, 1].set_title('D')

# Add legend to one of the subplots (they'll all have the same participants)
axs[0, 0].legend(handles=lines.values(), labels=lines.keys())

plt.tight_layout()
plt.show()


### Make a set of plots with regression lines
participants = sorted(participants)

model_params = pd.DataFrame(columns=['participant', 'intercept', 'slope'])

fig, axes = plt.subplots(2, 4, figsize=(16, 8))
axes = axes.flatten()

for index, participant in enumerate(participants):
    rows = []
    for label, series in latencies[participant]['difference'].items():
        for value in series:
            rows.append({'labels': label, 'differences': value})

    df = pd.DataFrame(rows)

    # Encode categories
    df['label_code'] = pd.Categorical(df['labels']).codes

    # Fit linear regression
    X = sm.add_constant(df['label_code'])  # adds intercept
    y = df['differences']
    model = sm.OLS(y, X).fit()

    # Predict values
    df['predicted'] = model.predict(X)

    # Extract slope and intercept
    intercept, slope = model.params

    model_params.loc[index] = {'participant': participant, 'intercept': intercept, 'slope': slope}

    ax = axes[index]
    ax.scatter(df['labels'], df['differences'])
    ax.plot(df['labels'], df['predicted'], color='red', label='Regression Line')

    ax.set_title(participant)
    ax.set_xlabel('Blocks')
    ax.set_ylabel('Difference [s]')

    # Add regression parameters as text in each subplot
    textstr = f"Slope: {slope:.4f}\nIntercept: {intercept:.4f}"
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
            fontsize=9, verticalalignment='top',
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.7))

plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

model_params.to_csv(Path(workingFolder, 'model_params_with_electrode_drift.csv'))
