import os
from ase import io
from ase.build import add_adsorbate
from ase.calculators.espresso import Espresso, EspressoProfile

default_pseudopotentials = {
    'Ti': 'Ti.pbe-spn-kjpaw_psl.1.0.0.UPF',
    'La': 'La.paw.z_11.atompaw.wentzcovitch.v1.2.upf',
    'Co': 'Co_pbe_v1.2.uspp.F.UPF',
    'O': 'O.pbe-n-kjpaw_psl.0.1.UPF',
    'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',
    'Sr': 'Sr_pbe_v1.uspp.F.UPF'
}

default_input_data = {
    'system': {
        'ecutwfc': 60,
        'ecutrho': 480,
        'occupations': 'smearing',
        'smearing': 'gaussian',
        'degauss': 0.07,
    },
    'control': {
        'calculation': 'scf',
        'restart_mode': 'from_scratch',
        'prefix': 'pwscf',
        'outdir': '/scratch',
        'disk_io': 'low',
        'verbosity': 'high'
    },
    'electrons': {
        'conv_thr': 1e-8,
        'mixing_mode': 'local-TF',
        'electron_maxstep': 200,
        'mixing_beta': 0.15,
        'diagonalization': 'david',
    }
}

adsorbates = {
    "OOH": [("O", 3), ("O", 3 + 1.245956), ("H", 3 + 1.245956 + .97907)],
    "OH": [("O", 3), ("H", 3 + .97907)],
    "O": [("O", 3)],
    "H": [("H", 3)],
}

class Generation:
    def __init__(self, cluster_name, cluster_command):
        self.cluster_name = cluster_name
        self.cluster_command = cluster_command

    def generate_files(self, material, pseudopotentials=None, input_data=None):
        if pseudopotentials is None:
            pseudopotentials = default_pseudopotentials
        if input_data is None:
            input_data = default_input_data

        input_data['control']['outdir'] = f'/bsuhome/{self.cluster_name}/scratch'

        profile = EspressoProfile(
            command=self.cluster_command,
            pseudo_dir=f"/bsuhome/{self.cluster_name}/q-e/pseudo"
        )

        header_folder = material[:-4]  # Create the header folder based on the material name without .xyz
        os.makedirs(header_folder, exist_ok=True)  # Ensure the header folder exists

        save_dirs = []

        read_material = io.read(material)

        for abs_name in adsorbates.keys():
            save_dir = os.path.join(header_folder, f"{material[:-4]}_{abs_name}")
            save_dirs.append(save_dir)
            for i in range(8):
                for j in range(8):
                    prev = read_material.copy()
                    for element, height in adsorbates[abs_name]:
                        add_adsorbate(read_material, element, height=height, position=(7.8254 * (i / 7), 7.8254 * (j / 7)))

                    bin_folder = "0" if i < 4 else "1"
                    dir_name = f"{i}-{j}"
                    file_name = f"{i}-{j}.xyz"
                    file_path = os.path.join(save_dir, bin_folder, dir_name)
                    os.makedirs(file_path, exist_ok=True)

                    io.write(os.path.join(file_path, file_name), read_material)

                    calc = Espresso(
                        profile=profile,
                        pseudopotentials=pseudopotentials,
                        input_data=input_data,
                        kpts=(4, 4, 1),
                        koffset=(0, 0, 0),
                    )
                    read_material.calc = calc
                    io.espresso.write_espresso_in(
                        os.path.join(file_path, "espresso.pwi"),
                        read_material,
                        input_data=input_data,
                        pseudopotentials=pseudopotentials,
                        kpts=(4, 4, 1),
                        koffset=(0, 0, 0)
                    )
                    read_material = prev.copy()

        return save_dirs