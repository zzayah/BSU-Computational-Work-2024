import os
import re
import pandas as pd

class Extract:
    def __init__(self, file_path, absorbates: list = None):
        self.file_path = file_path
        self.absorbates = absorbates if absorbates else ["O", "H", "OOH", "OH"]

        # Check if file_path is a filepath or just a filename
        if os.path.dirname(file_path):
            raise ValueError("Input file_path should be just a filename, not a filepath.")
        
        self.save_dirs = []
        # Extract the material name from the file path (e.g., SrCoO3 from SrCoO3.xyz)
        material_name = os.path.basename(file_path).replace('.xyz', '')
        for absorbate in self.absorbates:
            # Create directory names based on the material and absorbates
            self.save_dirs.append(f"{material_name}_{absorbate}")

    def extract_energy_from_file(self, file_path):
        with open(file_path, "r") as f:
            lines = f.readlines()
            for line in lines:
                if "!    total energy" in line:
                    # Extract the energy value using regex
                    energy = re.search(r"[-+]?\d*\.\d+|\d+", line)
                    if energy:
                        return float(energy.group(0))
        return None

    def process_folders(self, base_folder, start, end, df):
        for i in range(start, end):
            for j in range(start, end):
                folder_name = f"{i}-{j}"
                folder_path = os.path.join(base_folder, folder_name)
                if os.path.isdir(folder_path):
                    file_path = os.path.join(folder_path, "espresso.pwo")
                    if os.path.isfile(file_path):
                        energy = self.extract_energy_from_file(file_path)
                        if energy is not None:
                            df.loc[len(df)] = [i, j, energy]

    def extract_data(self):
        current_dir = os.getcwd()  # Get the current working directory
        for folder in self.save_dirs:
            # Initialize the dataframe
            df = pd.DataFrame(columns=["x", "y", "SCF"])

            # Process folders in "0"
            self.process_folders(os.path.join(folder, "0"), 0, 4, df)

            # Process folders in "1"
            self.process_folders(os.path.join(folder, "1"), 4, 8, df)

            # Save the dataframe to a CSV file in the current directory
            df.to_csv(os.path.join(current_dir, f"{folder}_OUT.csv"), index=False)
