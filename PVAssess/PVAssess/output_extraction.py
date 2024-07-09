import os
import re
import pandas as pd


def _extract_energy_from_file(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "!    total energy" in line:
                # Extract the energy value using regex
                energy = re.search(r"[-+]?\d*\.\d+|\d+", line)
                if energy:
                    return float(energy.group(0))
    return None


def _process_folders(base_folder, start, end, df):
    for i in range(start, end):
        for j in range(0, 8):
            folder_name = f"{i}-{j}"
            folder_path = os.path.join(base_folder, folder_name)
            if os.path.isdir(folder_path):
                file_path = os.path.join(folder_path, "espresso.pwo")
                if os.path.isfile(file_path):
                    energy = _extract_energy_from_file(file_path)
                    if energy is not None:
                        df.loc[len(df)] = [i, j, energy]


def Extract(file_path, absorbates=None):
    absorbates = absorbates if absorbates else ["O", "H", "OOH", "OH"]

    # Check if file_path is a filepath or just a filename
    if os.path.dirname(file_path):
        raise ValueError(
            "Inputted file_path should be the name of your material.")

    save_dirs = []
    # Extract the material name from the file path (e.g., SrCoO3 from SrCoO3.xyz)
    material_name = os.path.basename(file_path).replace('.xyz', '')
    for absorbate in absorbates:
        # Create directory names based on the material and absorbates
        save_dirs.append(f"{material_name}_{absorbate}")

    current_dir = os.getcwd()  # Get the current working directory
    for folder in save_dirs:
        # Initialize the dataframe
        df = pd.DataFrame(columns=["x", "y", "SCF"])

        # Process folders in "0"
        _process_folders(os.path.join(folder, "0"), 0, 4, df)

        # Process folders in "1"
        _process_folders(os.path.join(folder, "1"), 4, 8, df)

        # Save the dataframe to a CSV file in the current directory
        df.to_csv(os.path.join(current_dir, f"{folder}_OUT.csv"), index=False)
