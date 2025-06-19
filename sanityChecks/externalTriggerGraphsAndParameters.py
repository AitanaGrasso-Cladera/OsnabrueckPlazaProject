"""
This script is developed by Aziz Akkaya to create graphs and calculate parameters for the hardware events
"""
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
import pickle

# Settings
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 500)
pd.set_option('display.precision', 15)


# Define directories (delete the path)
workingFolder = Path('')
dataFolder = Path(workingFolder, 'Data')
triggerFolder = Path(workingFolder, 'triggerFiles')
eventFolder = Path(workingFolder, 'eventFiles')
hardwareFolder = Path(workingFolder, 'hardwareFiles')

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


# Load the pickle file
with open(Path(hardwareFolder, 'hardwareData.pkl'), 'rb') as f:
    latencies = pickle.load(f)

# External trigger names in the EEG data
target_labels = ['A_on', 'B_on', 'C_on', 'D_on']

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
plt.show(block=True)


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
