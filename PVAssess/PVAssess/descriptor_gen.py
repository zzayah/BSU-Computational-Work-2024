import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import itertools
import warnings
from abc import ABC
from typing import Callable, List, Union, Tuple
import ase as ase
from numpy.polynomial.chebyshev import Chebyshev as cheb
import numpy.typing as npt
from ase import Atoms
from ase.neighborlist import NeighborList
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


class CebychevCoefficients(ABC):
    order: int = 0
    radial: bool = True  # angular if False
    cutoff: float = 0.
    fc: Callable = None

    compositional: bool = False
    comp_weights: npt.NDArray[np.float64] = None
    structural: bool = False

    polynomia: List[Callable] = []

    def __init__(self, order: int = 1, compositional: bool = False, structural: bool = False, radial: bool = True, cutoff: float = 0.) -> None:
        self.order = order
        self.compositional = compositional
        self.structural = structural
        self.radial = radial
        if self.radial:
            if cutoff == 0:
                raise ValueError()
            self.cutoff = cutoff
        else:
            self.cutoff = np.pi
        self._generate_functions()

    def _generate_functions(self) -> None:
        self.fc = wraprcut(cosfc, self.cutoff)
        for i in range(self.order):
            pol_i = wraprcut(wraporder(dual_basis_function, i), self.cutoff)
            self.polynomia.append(pol_i)

    def _init_weights(self, types: npt.NDArray[np.int64]) -> None:
        # This will have to be implemented later with a correct weight function
        return NotImplementedError()

    def calculate(self, distances: npt.NDArray, vectors: npt.NDArray = None, types: npt.NDArray[np.int64] = None):
        if self.compositional:
            self._init_weights(types)
        if (self.radial):
            result = self._calculate_radial(distances)
            return result
        else:
            result = self._calculate_angular(distances, vectors)
            return result

    def _calculate_radial(self, distances):
        fci = self.fc(distances)
        if self.structural:
            structural_coefficients = np.zeros(self.order)
        if self.compositional:
            composition_coefficients = np.zeros(self.order)

        for j in fci.nonzero()[0]:
            for a in range(self.order):
                rij = distances[j]
                if self.structural:
                    structural_coefficients[a] += self.polynomia[a](
                        rij) * fci[j]
                if self.compositional:
                    composition_coefficients[a] += self.polynomia[a](
                        rij) * fci[j] * self.comp_weights[a]

        if self.structural:
            if self.compositional:
                return structural_coefficients, composition_coefficients
            return structural_coefficients
        if self.compositional:
            return composition_coefficients

    def _calculate_angular(self, distances, vectors):
        fci = self.fc(distances)
        if self.structural:
            structural_coefficients = np.zeros(self.order)
        if self.compositional:
            composition_coefficients = np.zeros(self.order)

        # iterate through distance pairs and their corresponding vector pairs
        for j in fci.nonzero()[0]:
            rij = distances[j]
            rij_vec = vectors[j]
            for k in fci.nonzero()[0]:
                if k < j:
                    continue
                rik = distances[k]
                rik_vec = vectors[k]
                angle = np.arccos(
                    np.clip(np.dot(rij_vec, rik_vec) / (rij * rik), -1, 1))
                for a in range(self.order):
                    if self.structural:
                        structural_coefficients[a] += self.polynomia[a](
                            angle) * fci[j] * fci[k]
                    if self.compositional:
                        composition_coefficients[a] += self.polynomia[a](
                            angle) * fci[j] * fci[k] * self.comp_weights[a]

        if self.structural:
            if self.compositional:
                return structural_coefficients, composition_coefficients
            return structural_coefficients
        if self.compositional:
            return composition_coefficients

    @staticmethod
    def _generate_vectors_and_distances(atoms) -> Tuple[npt.NDArray, npt.NDArray]:
        positions = atoms.get_positions()
        n_atoms = len(positions)
        distances = np.zeros((n_atoms, n_atoms))
        vectors = np.zeros((n_atoms, n_atoms, 3))

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                vec = positions[j] - positions[i]
                dist = np.linalg.norm(vec)
                distances[i, j] = dist
                distances[j, i] = dist
                vectors[i, j] = vec
                vectors[j, i] = -vec

        return distances, vectors


def wraprcut(
    f: Callable,
    rcut: float,
) -> Callable:
    def wrapped(*args) -> float:
        return f(rcut, *args)

    return wrapped


def wraporder(
    f: Callable,
    order: int,
) -> Callable:
    def wrapped(*args) -> float:
        return f(order, *args)

    return wrapped


# def wrapatoms(
#     f: Callable,
#     atoms: ase.Atoms,
# ) -> Callable:
#     def wrapped(*args) -> float:
#         return f(atoms, *args)

#     return wrapped


def dual_basis_function(
        order: int,
        rcut: float,
        x: float):
    return cheb.basis(order)(2 * x/rcut - 1)


def Descriptors():

    from ase.io import read
    # the coordinates in the file are in units of alat
    mos2_defect = read('SrCoO3_A.xyz')
    # alat is in a.u. (bohr), need to convert to Angstrom
    alat = 17.8707*0.5291772
    mos2_defect.positions = mos2_defect.positions * \
        alat  # convert positions to Angstrom
    mos2_defect.cell = np.array(
        [1., 1.155671, 2.643601])*alat  # set system cell
    mos2_defect.pbc = [1, 1, 1]  # and add boundary conditions

    elements = mos2_defect.get_chemical_symbols()
    types = sorted(list(dict.fromkeys(elements)))
    itypes = np.array([types.index(e) for e in elements]).astype('int')
    icombos = list(itertools.combinations_with_replacement(
        np.unique(itypes), 2))
    # print(len(elements))
    # print(itypes)
    # print(icombos)

    relax_data = pd.read_csv('SrCoO3_A_H_OUT.csv', index_col=0)
    # we may want to check if using all the data is ok
    filtered_data = relax_data[relax_data['SCF'] < -4399].copy()
    # filtered_data=filtered_data.drop(columns=['Fx','Fy','Fz'])
    # print(filtered_data)

    filtered_data = filtered_data.reset_index()
    virtual_atoms = np.array(np.column_stack(
        (filtered_data['x'], filtered_data['y'], filtered_data['z'])))
    virtual_atoms_indexes = filtered_data.index
    n_virtual_atoms = virtual_atoms.shape[0]

    rcuts = [5.]
    order = 15

    # RDF w/ Chebyshev coefficients
    total_atoms = mos2_defect.copy()
    total_atoms.append('He')
    for i in range(n_virtual_atoms):
        total_atoms.positions[-1] = virtual_atoms[i]
        distances = total_atoms.get_distances(
            len(total_atoms)-1, list(range(len(total_atoms)-1)), mic=True)
        for rcut in rcuts:
            cheb = CebychevCoefficients(
                order=order, compositional=False, structural=True, radial=True, cutoff=rcut)
            coefficients = cheb.calculate(distances=distances)
            for a in range(order):
                filtered_data.loc[virtual_atoms_indexes[i],
                                  f'rdf_rcut_{rcut}_order_{a}'] = coefficients[a]

    filtered_data.to_csv('Chebyshev_RDF_MoS2.csv', index=False)

    # # Extract RDF data for the first rcut and order
    # rdf_column = f'rdf_rcut_{rcuts[0]}_order_{0}'  # Change 0 to the desired order if necessary
    # rdf_data = filtered_data[rdf_column]

    # # Plot the histogram
    # plt.hist(rdf_data, bins=50, edgecolor='black', alpha=0.9, color='orange')
    # plt.xlabel('RDF Value', labelpad=10, fontsize="25")
    # plt.ylabel('Frequency', labelpad=10, fontsize="25")
    # #plt.title(f'Chebyshev Approximation for RDF', fontsize="20")  # Change 0 to the desired order if necessary
    # plt.xticks([])
    # plt.yticks([])
    # plt.show()

    def compute_rdf(distances, rcut, bin_width=0.1):
        bins = np.arange(0, rcut + bin_width, bin_width)
        rdf, _ = np.histogram(distances, bins=bins)
        return rdf

    total_atoms = mos2_defect.copy()
    total_atoms.append('He')

    for i in range(n_virtual_atoms):
        total_atoms.positions[-1] = virtual_atoms[i]

        # Create NeighborList to find distances
        cutoffs = [max(rcuts)] * len(total_atoms)
        neighbor_list = NeighborList(
            cutoffs, skin=0.0, self_interaction=False, bothways=True)
        neighbor_list.update(total_atoms)

        indices, offsets = neighbor_list.get_neighbors(len(total_atoms) - 1)
        distances = total_atoms.get_distances(
            len(total_atoms) - 1, indices, mic=True)

        for rcut in rcuts:
            rdf = compute_rdf(distances, rcut)

            for a, value in enumerate(rdf):
                filtered_data.loc[virtual_atoms_indexes[i],
                                  f'rdf_rcut_{rcut}_bin_{a}'] = value

    # # Extract RDF data for the first rcut and order
    # rdf_column = f'rdf_rcut_{rcuts[0]}_order_{0}'  # Change 0 to the desired order if necessary
    # rdf_data = filtered_data[rdf_column]

    # # Plot the histogram
    # plt.hist(rdf_data, bins=50, edgecolor='black', alpha=0.9, color='blue')
    # plt.xlabel('RDF Value', labelpad=10, fontsize="25")
    # plt.ylabel('Frequency', labelpad=10, fontsize="25")
    # # plt.title(f'Calculated RDF', fontsize="20")  # Change 0 to the desired order if necessary
    # plt.xticks([])
    # plt.yticks([])
    # plt.show()

    total_atoms = mos2_defect.copy()
    total_atoms.append('He')

    for i in range(n_virtual_atoms):
        total_atoms.positions[-1] = virtual_atoms[i]
        distances = total_atoms.get_distances(
            len(total_atoms)-1, list(range(len(total_atoms)-1)), mic=True)
        vectors = total_atoms.get_distances(
            len(total_atoms)-1, list(range(len(total_atoms)-1)), mic=True, vector=True)
        for rcut in rcuts:
            cheb = CebychevCoefficients(
                order=order, compositional=False, structural=True, radial=False, cutoff=rcut)
            coefficients = cheb.calculate(distances=distances, vectors=vectors)
            for a in range(order):
                filtered_data.loc[virtual_atoms_indexes[i],
                                  f'adf_rcut_{rcut}_order_{a}'] = coefficients[a]

    # for i in [0]:
    #     total_atoms.positions[-1] = virtual_atoms[i]
    #     distances = total_atoms.get_distances(len(total_atoms)-1, list(range(len(total_atoms)-1)), mic=True)
    #     plt.hist(distances,bins=100)

    return None
