import os
from ase import io
from ase.build import add_adsorbate
from ase.calculators.espresso import Espresso, EspressoProfile

def Generate(material, cluster_name="Borah", cluster_username=None, cluster_command=None, psuedo_dir=None, scratch_dir=None, pseudopotentials=None, input_data=None):

    if cluster_username is None:
        raise ValueError("Please specify your cluster username.")
    if cluster_command is None:
        print("As of 7/2/24, the cluster_command perameter for the Generate function may be == "", though this is subject to change in future QE or ASE versions. If you choose to leave this field empty, please ensure that the cluster command can be seen in the .pwi files.")
        # raise ValueError("Please specify your cluster command.")
    if psuedo_dir is None and cluster_name == "Borah":
        psuedo_dir = f"/bsuhome/{cluster_username}/q-e/pseudo"
    elif psuedo_dir is None:
        raise ValueError("Please specify the path to the pseudopotentials directory in your cluster directory.")
    if scratch_dir is None and cluster_name == "Borah":
        scratch_dir = f"/bsuhome/{cluster_username}/scratch"
    elif scratch_dir is None:
        raise ValueError("Please specify the path to the scratch directory in your cluster directory.")

    # All SSSP_acc_PBE Psuedopotentials for QE
    default_pseudopotentials = {
        'Ag': 'ag_pbe_v1.4.uspp.F.UPF',
        'Al': 'Al.pbe-n-kjpaw_psl.1.0.0.UPF',
        'Ar': 'Ar.pbe-n-rrkjus_psl.1.0.0.UPF',
        'As': 'As.pbe-n-rrkjus_psl.0.2.UPF',
        'Au': 'Au_ONCV_PBE-1.0.upf',
        'B': 'B.pbe-n-kjpaw_psl.0.1.UPF',
        'Ba': 'Ba_ONCV_PBE-1.0.upf',
        'Be': 'Be_ONCV_PBE-1.0.upf',
        'Bi': 'Bi.pbe-dn-kjpaw_psl.0.2.2.UPF',
        'Br': 'br_pbe_v1.4.uspp.F.UPF',
        'Ca': 'Ca_pbe_v1.uspp.F.UPF',
        'Cd': 'Cd.pbe-dn-rrkjus_psl.0.3.1.UPF',
        'Ce': 'Ce.GGA-PBE-paw-v1.0.UPF',
        'Cl': 'Cl.pbe-n-rrkjus_psl.1.0.0.UPF',
        'Co': 'Co_pbe_v1.2.uspp.F.UPF',
        'Cr': 'cr_pbe_v1.5.uspp.F.UPF',
        'Cs': 'Cs_pbe_v1.uspp.F.UPF',
        'Cu': 'Cu_pbe_v1.2.uspp.F.UPF',
        'Dy': 'Dy.GGA-PBE-paw-v1.0.UPF',
        'Er': 'Er.GGA-PBE-paw-v1.0.UPF',
        'Eu': 'Eu.GGA-PBE-paw-v1.0.UPF',
        'Fe': 'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF',
        'F': 'f_pbe_v1.4.uspp.F.UPF',
        'Ga': 'Ga.pbe-dn-kjpaw_psl.1.0.0.UPF',
        'Gd': 'Gd.GGA-PBE-paw-v1.0.UPF',
        'Ge': 'Ge.pbe-dn-kjpaw_psl.1.0.0.UPF',
        'H': 'H.pbe-rrkjus_psl.1.0.0.UPF',
        'Hf': 'Hf.pbe-spdfn-kjpaw_psl.1.0.0.UPF',
        'Hg': 'Hg_pbe_v1.uspp.F.UPF',
        'Ho': 'Ho.GGA-PBE-paw-v1.0.UPF',
        'In': 'In.pbe-dn-rrkjus_psl.0.2.2.UPF',
        'Ir': 'Ir_pbe_v1.2.uspp.F.UPF',
        'I': 'I_pbe_v1.uspp.F.UPF',
        'K': 'K.pbe-spn-rrkjus_psl.1.0.0.UPF',
        'Kr': 'Kr.pbe-n-rrkjus_psl.0.2.3.UPF',
        'La': 'La.paw.z_11.atompaw.wentzcovitch.v1.2.upf',
        'Li': 'li_pbe_v1.4.uspp.F.UPF',
        'Lu': 'Lu.GGA-PBE-paw-v1.0.UPF',
        'Mg': 'mg_pbe_v1.4.uspp.F.UPF',
        'Mn': 'Mn.pbe-spn-kjpaw_psl.0.3.1.UPF',
        'Mo': 'Mo_ONCV_PBE-1.0.upf',
        'Na': 'Na_pbe_v1.uspp.F.UPF',
        'Nb': 'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF',
        'Nd': 'Nd.GGA-PBE-paw-v1.0.UPF',
        'Ne': 'Ne.pbe-n-kjpaw_psl.1.0.0.UPF',
        'Ni': 'ni_pbe_v1.4.uspp.F.UPF',
        'O': 'O.pbe-n-kjpaw_psl.0.1.UPF',
        'Os': 'Os.pbe-spfn-rrkjus_psl.1.0.0.UPF',
        'P': 'P.pbe-n-rrkjus_psl.1.0.0.UPF',
        'Pb': 'Pb.pbe-dn-kjpaw_psl.0.2.2.UPF',
        'Pd': 'Pd.pbe-spn-kjpaw_psl.1.0.0.UPF',
        'Pm': 'Pm.GGA-PBE-paw-v1.0.UPF',
        'Po': 'Po.pbe-dn-rrkjus_psl.1.0.0.UPF',
        'Pr': 'Pr.GGA-PBE-paw-v1.0.UPF',
        'Pt': 'Pt.pbe-spfn-rrkjus_psl.1.0.0.UPF',
        'Rb': 'Rb_ONCV_PBE-1.0.upf',
        'Re': 'Re_pbe_v1.2.uspp.F.UPF',
        'Rh': 'Rh.pbe-spn-kjpaw_psl.1.0.0.UPF',
        'Rn': 'Rn.pbe-dn-rrkjus_psl.1.0.0.UPF',
        'Ru': 'Ru_ONCV_PBE-1.0.upf',
        'Sb': 'sb_pbe_v1.4.uspp.F.UPF',
        'Sc': 'Sc_pbe_v1.uspp.F.UPF',
        'Se': 'Se_pbe_v1.uspp.F.UPF',
        'Si': 'Si.pbe-n-rrkjus_psl.1.0.0.UPF',
        'Sm': 'Sm.GGA-PBE-paw-v1.0.UPF',
        'Sn': 'Sn_pbe_v1.uspp.F.UPF',
        'Sr': 'Sr_pbe_v1.uspp.F.UPF',
        'S': 'S_pbe_v1.2.uspp.F.UPF',
        'Ta': 'Ta.pbe-spfn-rrkjus_psl.1.0.0.UPF',
        'Tb': 'Tb.GGA-PBE-paw-v1.0.UPF',
        'Tc': 'Tc_ONCV_PBE-1.0.upf',
        'Te': 'Te_pbe_v1.uspp.F.UPF',
        'Ti': 'Ti.pbe-spn-kjpaw_psl.1.0.0.UPF',
        'Tl': 'Tl.pbe-dn-rrkjus_psl.1.0.0.UPF',
        'Tm': 'Tm.GGA-PBE-paw-v1.0.UPF',
        'V': 'V_pbe_v1.uspp.F.UPF',
        'W': 'W_pbe_v1.2.uspp.F.UPF',
        'Xe': 'Xe.pbe-dn-rrkjus_psl.1.0.0.UPF',
        'Yb': 'Yb.GGA-PBE-paw-v1.0.UPF',
        'Y': 'Y_pbe_v1.uspp.F.UPF',
        'Zn': 'Zn_pbe_v1.uspp.F.UPF',
        'Zr': 'Zr_pbe_v1.uspp.F.UPF'
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

    if pseudopotentials is None:
        pseudopotentials = default_pseudopotentials
    if input_data is None:
        input_data = default_input_data

    input_data['control']['outdir'] = f'/bsuhome/{cluster_username}/scratch'

    profile = EspressoProfile(
        command=cluster_command,
        pseudo_dir=psuedo_dir,
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