import numpy as np
import ase
from ase import Atoms
from ase.build import make_supercell
from ase.visualize import view
from ase.io import read, write
from ase import io

def ABSite(action, filepath):
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
        AtomsLayerB_R = AllLayers[0] + "_SB"
    elif VisualizeLayers[len(VisualizeLayers) - 1] == 'B':
        AtomsLayerB_R = AllLayers[len(VisualizeLayers) - 1] + "_SB"
    
    if VisualizeLayers[0] == 'A':
        AtomsLayerA_R = AllLayers[0] + "_SB"
    if VisualizeLayers[len(VisualizeLayers) - 1] == 'A':
        AtomsLayerA_R = AllLayers[len(VisualizeLayers) - 1] + "_SB"
    


    NumberOfLayers = len(AllLayers)
    print(AllLayers)
    print(VisualizeLayers)
    #it has the perfect layers to delete
    AtomNumber = []
    for layer in AllLayers:
        globals()[layer+"_SB"] = []
        
    for atom in read(SiteBase):
        ZName = "Z_" + str(atom.position[2]) + "_SB"
        globals()[ZName].append([atom.symbol, atom.position[0], atom.position[1], atom.index])
    
    for layer in AllLayers:
        LName = layer+"_SB"
        print(LName, globals()[LName])    
    NumberOfAtoms = len(globals()[layer+"_SB"])
    AtomNumber.append(NumberOfAtoms)
    print(AtomNumber)

    FileA = "ASite_"+FileName+".xyz"
    FileB = "BSite_"+FileName+".xyz"
    def CreateA():
        # return FileName, AtomNumber, NumberOfLayers
        with open(FileA, 'w') as f:
            f.write('')
        AtomsOnLayer = globals()[AtomsLayerB_R]
        atomNumbersL = len(read(SiteBase))-len(AtomsOnLayer)

        with open('General_' + FileName + '.xyz', 'r') as f:
            for line in f:
                if 'Lattice=' in line:
                    lattice_info = line.strip()
        with open(FileA, 'a') as f:
            f.write(f"{atomNumbersL}\n")
        with open(FileA, 'a') as f:
            f.write(f"{lattice_info}\n")

        for atom in ViewBase:
            if [atom.symbol, atom.position[0], atom.position[1], atom.index] not in AtomsOnLayer:
                with open(FileA, 'a') as f:
                    f.write(f"{atom_to_xyz_line(atom)}")
                # Read the atomic structure from XYZ file
        atomsA = read(FileA)

        # Get the current cell parameters
        cell_paramsA = atomsA.get_cell()

        # Modify the cell dimensions to change vacuum size (adjust along the z-axis)
        vacuum_size = 10.0  # Adjust this value to change the size of the vacuum
        cell_paramsA[2, 2] += vacuum_size  # Increase or decrease vacuum size along z-axis

        # Calculate current middle along the z-axis
        middle_zA = (cell_paramsA[2, 2] + cell_paramsA[0, 2]) / 2.0

        # Calculate shift needed to center the structure
        shift_zA = middle_zA - np.mean(atomsA.positions[:, 2])


        # Shift positions of all atoms to center the vacuum
        atomsA.positions[:, 2] += shift_zA

        # Set the modified cell parameters back to the structure
        atomsA.set_cell(cell_paramsA)

        # Visualize the structure with modified vacuum size and centered along the z-axis
        view(atomsA)
   
    
    def CreateB():
        with open(FileB, 'w') as f:
            f.write('')
        AtomsOnLayer = globals()[AtomsLayerA_R]
        atomNumbersL = len(read(SiteBase))-len(AtomsOnLayer)

        with open('General_' + FileName + '.xyz', 'r') as f:
            for line in f:
                if 'Lattice=' in line:
                    lattice_info = line.strip()
        with open(FileB, 'a') as f:
            f.write(f"{atomNumbersL}\n")
        with open(FileB, 'a') as f:
            f.write(f"{lattice_info}\n")

        for atom in ViewBase:
            if [atom.symbol, atom.position[0], atom.position[1], atom.index] not in AtomsOnLayer:
                with open(FileB, 'a') as f:
                    f.write(f"{atom_to_xyz_line(atom)}")
            # Read the atomic structure from XYZ file
        
        atomsB = read(FileB)
        # Get the current cell parameters
    
        cell_paramsB = atomsB.get_cell()
        # Modify the cell dimensions to change vacuum size (adjust along the z-axis)
        vacuum_size = 10.0  # Adjust this value to change the size of the vacuum
  
        cell_paramsB[2, 2] += vacuum_size
        # Calculate current middle along the z-axis
     
        middle_zB = (cell_paramsB[2, 2] + cell_paramsB[0, 2]) / 2.0
        # Calculate shift needed to center the structure
 
        shift_zB = middle_zB - np.mean(atomsB.positions[:, 2])

        # Shift positions of all atoms to center the vacuum

        atomsB.positions[:, 2] += shift_zB
        # Set the modified cell parameters back to the structure

        atomsB.set_cell(cell_paramsB)
        # Visualize the structure with modified vacuum size and centered along the z-axis

        view(atomsB)
    
    
    
    if action == "Both":
        CreateA()
        CreateB()
    if action == "A":
        CreateA()
    if action == "B":
        CreateA()




def Create_Structs(filepath, write_to: str = "./"):
    action = "Both"
    ABSite(action, filepath)
    return # returns nothing. writes A & B site .xyz strutures to specifced "write_to". If write_to not defied, writes to current diretory (./).


def Create_Struct_A(filepath, write_to: str = "./"):
    action = "A"
    ABSite(action, filepath)
    return # returns nothing. writes A site .xyz struture to specifced "write_to". If write_to not defied, writes to current diretory (./).


def Create_Struct_B(filepath, write_to: str = "./"):
    action = "B"
    ABSite(action, filepath)
    return # returns nothing. writes B site .xyz struture to specifced "write_to". If write_to not defied, writes to current diretory (./).


