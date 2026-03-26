### checks the sampling rate of the scene camera

import os
import pandas as pd
import numpy as np

base_path = ""

sampling_rates = []

for participant_folder in os.listdir(base_path):
    if participant_folder.startswith("Participant"):
        participant_number = participant_folder.split()[-1]
        participant_id = f"P{participant_number}"
        folder_path = os.path.join(base_path, participant_folder)

        for recording_folder in os.listdir(folder_path):
            if recording_folder.startswith(f"P{participant_number}_"):
                recording_number = recording_folder.split('_')[-1]
                recording_path = os.path.join(folder_path, recording_folder)

                try:
                    files = {
                        "world_timestamps": os.path.join(recording_path, "world_timestamps.csv"),
                    }
                    
                    # Load data
                    data = {name: pd.read_csv(file) for name, file in files.items()}
                    df_ts = data["world_timestamps"]

                    # sampling rate
                    duration_ns = df_ts['timestamp [ns]'].iloc[-1] - df_ts['timestamp [ns]'].iloc[0]
                    duration_sec = duration_ns / 1e9
                    
                    if duration_sec > 0:
                        sampling_rate = len(df_ts) / duration_sec
                        sampling_rates.append(sampling_rate)
                        
                except Exception as e:
                    print(f"Skipping {recording_folder}: {e}")

print(sampling_rates)

if sampling_rates:
    stats = {
        "Mean": np.mean(sampling_rates),
        "SD":   np.std(sampling_rates),
        "Min":  np.min(sampling_rates),
        "Max":  np.max(sampling_rates)
    }

    for metric, value in stats.items():
        print(f"{metric}: {value:.4f} Hz")
