import numpy as np
import ase
from ase import Atoms
from ase.build import make_supercell
from ase.visualize import view
from ase.io import read, write
from ase import io

def process_cif_file(filepath):
    # Extract filename without extension
    FileName = filepath.split('.')[0]
    FileType = filepath.split('.')[-1]

    File = read(filepath)
    
    CompoundA = ''
    CompoundB = ''
    Move = 0
    for letter in FileName:
        if letter.isupper():
            Move += 1
        if Move == 1:
            CompoundA = CompoundA + letter
        if Move == 2:
            CompoundB = CompoundB + letter

    write(FileName + ".xyz", File)
    File = FileName + ".xyz"

    atom = read(File)
    rotation_offset = 0
    rotate_scale = 90

    supercell_atoms = make_supercell(atom, [[3, 0, 0], [0, 3, 0], [0, 0, 3]])
    corner_atoms = []

    while CompoundA not in corner_atoms:
        corner_atoms = []
        supercell_atoms.rotate(rotation_offset, 'z')
        
        min_x = np.max(supercell_atoms.positions[:, 0])
        min_y = np.max(supercell_atoms.positions[:, 1])
        max_z = np.min(supercell_atoms.positions[:, 2])

        for atom in supercell_atoms:
            if atom.position[0] <= min_x:
                if f"{atom.position[0]:.8f}" == "-0.00000000":
                    min_x = -0.00000000
                else:
                    min_x = atom.position[0]
            if atom.position[1] <= min_y:
                if f"{atom.position[1]:.8f}" == "-0.00000000":
                    min_y = -0.00000000
                else:
                    min_y = atom.position[1]

        for atom in supercell_atoms:
            if f"{atom.position[0]:.8f}" == f"{min_x:.8f}" and f"{atom.position[1]:.8f}" == f"{min_y:.8f}":
                corner_atoms.append(atom.symbol)
        
        if rotation_offset > 0:
            displacement = [-min_x, -min_y, 0]
            supercell_atoms.translate(displacement)
        rotation_offset = rotation_offset + rotate_scale

    write('SS_' + FileName + '.xyz', supercell_atoms)

    Bush = read('SS_' + FileName + '.xyz')
    lattice_info = None
    General = 'General_' + FileName + '.xyz'
    with open(General, 'w') as f:
        f.write('')

    with open('SS_' + FileName + '.xyz', 'r') as f:
        for line in f:
            if 'Lattice=' in line:
                lattice_info = line.strip()

    AllElements = [atom.symbol for atom in Bush]
    TypesOfElements = []
    for Elements in AllElements:
        if Elements not in TypesOfElements:
            TypesOfElements.append(Elements)
    ElementCorrespondentList = {}
    for Element in TypesOfElements:
        with open('SS_' + FileName + '.xyz', 'r') as f:
            for line in f:
                if Element in line:
                    Number = line[len(line) - 2]
                    ElementCorrespondentList[Element] = Number

    def atom_to_xyz_line(atom):
        symbol = atom.symbol
        coordinates = atom.position
        if len(symbol) == 2:
            line = f"{symbol}      {coordinates[0]:.8f}         {coordinates[1]:.8f}        {coordinates[2]:.8f}         {ElementCorrespondentList[atom.symbol]}\n"
        else:
            line = f"{symbol}       {coordinates[0]:.8f}         {coordinates[1]:.8f}        {coordinates[2]:.8f}         {ElementCorrespondentList[atom.symbol]}\n"
        return line

    atomNumbers = 0
    for atom in Bush:
        x = atom.position[0]
        y = atom.position[1]
        z = atom.position[2]
        if x != 0 and y != 0:
            atomNumbers = atomNumbers + 1

    with open(General, 'a') as f:
        f.write(f"{atomNumbers}\n")
    with open(General, 'a') as f:
        f.write(f"{lattice_info}\n")

    for atom in Bush:
        x = atom.position[0]
        y = atom.position[1]
        z = atom.position[2]
        if x != 0 and y != 0:
            with open(General, 'a') as f:
                xyz_line = atom_to_xyz_line(atom)
                f.write(xyz_line)

    GeneralFile = read(General)
    positions = GeneralFile.get_positions()
    min_coords = np.min(GeneralFile.positions, axis=0)
    max_coords = np.max(GeneralFile.positions, axis=0)

    mol_size_x = max_coords[0] - min_coords[0]
    mol_size_y = max_coords[1] - min_coords[1]
    mol_size_z = max_coords[2] - min_coords[2]

    vacuum_size_x = max_coords[0] - min_coords[0]
    vacuum_size_y = max_coords[1] - min_coords[1]
    vacuum_size_z = max_coords[2] - min_coords[2]

    center_x_rotated = min_coords[0]
    center_y_rotated = min_coords[1]
    center_z_rotated = min_coords[2]
    GeneralFile.positions -= [center_x_rotated, center_y_rotated, center_z_rotated]

    GeneralFile.set_cell([vacuum_size_x, vacuum_size_y, vacuum_size_z])
    write(General, GeneralFile)

    TypeOfZ = []
    for atom in GeneralFile:
        ULI = atom.position[2]
        if ULI not in TypeOfZ:
            TypeOfZ.append(ULI)

    for Z in TypeOfZ:
        ZLayerName = "Z_" + str(Z)
        globals()[ZLayerName] = []
        for atom in GeneralFile:
            if atom.position[2] == Z:
                if [atom.symbol, atom.position[0], atom.position[1], atom.index] not in globals()[ZLayerName]:
                    globals()[ZLayerName].append([atom.symbol, atom.position[0], atom.position[1], atom.index])

    GreatestXVal = 0
    GreatestYVal = 0
    IndexesToDelete = []
    for atom in GeneralFile:
        if atom.position[0] > GreatestXVal:
            GreatestXVal = atom.position[0]
        if atom.position[1] > GreatestXVal:
            GreatestYVal = atom.position[1]

    for Z in TypeOfZ:
        ZLayerName = "Z_" + str(Z)
        for atom_info in globals()[ZLayerName]:
            if GreatestXVal in atom_info and atom_info[3] not in IndexesToDelete:
                IndexesToDelete.append(atom_info[3])

    atomNumbersG = 0
    for atom in GeneralFile:
        atomNumbersG = atomNumbersG + 1
    SiteBase = "SB_" + FileName + ".xyz"
    with open(SiteBase, 'w') as f:
        f.write('')
    atomNumbersG = atomNumbersG - len(IndexesToDelete)

    with open('General_' + FileName + '.xyz', 'r') as f:
        for line in f:
            if 'Lattice=' in line:
                lattice_info = line.strip()
    with open(SiteBase, 'a') as f:
        f.write(f"{atomNumbersG}\n")
    with open(SiteBase, 'a') as f:
        f.write(f"{lattice_info}\n")

    for atom in GeneralFile:
        if atom.index not in IndexesToDelete:
            with open(SiteBase, 'a') as f:
                f.write(f"{atom_to_xyz_line(atom)}")

    ViewBase = read(SiteBase)

    AllLayers = []
    for Z in TypeOfZ:
        ZLayerName = "Z_" + str(Z)
        AllLayers.append(ZLayerName)
    AllLayers = sorted(AllLayers)

    BSiteLayers = []
    for Z in TypeOfZ:
        ZLayerName = "Z_" + str(Z)
        for atominfo in globals()[ZLayerName]:
            if atominfo[0] == CompoundB and ZLayerName not in BSiteLayers:
                BSiteLayers.append(ZLayerName)

    ASiteLayers = []
    for Z in TypeOfZ:
        ZLayerName = "Z_" + str(Z)
        for atominfo in globals()[ZLayerName]:
            if atominfo[0] == CompoundA and ZLayerName not in ASiteLayers:
                ASiteLayers.append(ZLayerName)
    VisualizeLayers = []
    for Layer in AllLayers:
        if Layer in BSiteLayers:
            VisualizeLayers.append("B")
        elif Layer in ASiteLayers:
            VisualizeLayers.append("A")

    if VisualizeLayers[0] == 'B':
        AtomsLayer = AllLayers[0] + "_SB"
    elif VisualizeLayers[len(VisualizeLayers) - 1] == 'B':
        AtomsLayer = AllLayers[len(VisualizeLayers) - 1] + "_SB"
    for i in range(0, len(VisualizeLayers)):
        if VisualizeLayers[i] == 'B':
            if VisualizeLayers[i - 1] == 'A' and VisualizeLayers[i + 1] == 'A':
                AtomsLayer = AllLayers[i] + "_SB"

    NumberOfLayers = len(AllLayers)

    AtomNumber = []
    for layer in AllLayers:
        NumberOfAtoms = len(globals()[layer])
        AtomNumber.append(NumberOfAtoms)

    # return FileName, AtomNumber, NumberOfLayers

# Example usage:
filepath = 'SrTiO3.cif'
process_cif_file(filepath)

