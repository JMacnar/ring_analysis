import sys, os, requests
import numpy as np
from functools import reduce
import argparse
import gzip


sys.path.append(os.path.join(os.path.dirname(__file__), 'bin'))
os.environ["BIOSHELL_DATA_DIR"] = os.path.join(os.path.dirname(__file__), 'data')

from pybioshell.core.data.io import Pdb, Cif
from pybioshell.core.chemical import PdbMolecule
from pybioshell.core.chemical import find_rings
from pybioshell.core.data.structural.selectors import SelectResidueByName, ChainSelector, SelectChainResidues
from pybioshell.core.calc.structural import SaturatedRing6Geometry
from pybioshell.core.data.structural import Residue, Chain, Structure
from pybioshell.core.chemical import MonomerStructure

#from pybioshell.core import BioShellVersion
#print(BioShellVersion().to_string())


def extract_ligand(pdb_file_name, ligand_name, cutoff_distance):
    #check if the file is not empty
    if os.stat(pdb_file_name).st_size == 0:
        print(pdb_file_name, "has", os.stat(pdb_file_name).st_size)
        return
    #read file as bioshell structure without hydrogens and alternative atom positions
    if os.stat(pdb_file_name).st_size == 0:
        print(pdb_file_name, "has", os.stat(pdb_file_name).st_size)
        return
    s = Pdb(str(pdb_file_name), "is_not_hydrogen is_not_alternative", False, False).create_structure(0)
    ligands = []
    for ic in range(s.count_chains()):
        chain = s[ic]
        for ir in range(chain.count_residues()):
            if chain[ir].residue_type().code3 == ligand_name:
                ligands.append(chain[ir])
    for l in ligands:
        fname = "%s-%d-%s-%s.pdb" % (l.residue_type().code3, l.id(), l.owner().id(), s.code())
        fout = open(fname, "w")

        for ic in range(s.count_chains()):
            chain = s[ic]
            for ir in range(chain.count_residues()):
                r = chain[ir]
                if r.min_distance(l) < cutoff_distance:
                    if r.id() == l.id():
                        for ai in range(r.count_atoms()):
                            fout.write(r[ai].to_pdb_line() + "\n")
                    if r.residue_type().code3 != l.residue_type().code3:
                        for ai in range(r.count_atoms()):
                            fout.write(r[ai].to_pdb_line() + "\n")

        #full CONNECT section from pdb file is added to the ligand's pdb - it is crucial to define rings
        for line in open(pdb_file_name).readlines():
            if line.startswith("CONECT"):
                fout.write(line)

        fout.close()
        return l, fname


def download_ideal_cif(ligand_code):
    try:
        r = requests.get("https://files.rcsb.org/ligands/download/%s.cif" % ligand_code, allow_redirects=True)
        #r = requests.get("https://files.rcsb.org/ligands/download/%s_ideal.sdf" % ligand_code, allow_redirects=True)
        open(ligand_code + ".cif", "wb").write(r.content)
    except Exception as error:
        print("An issue occurred during", ligand_code, "downloading", error)
    mm = MonomerStructure.from_cif(ligand_code + ".cif")
    return mm

def load_atoms_counts(fname):
    f = open(fname)
    for i in range(3): f.readline()
    tokens = f.readline().strip().split()
    n = int(tokens[0][0:3])
    number_of_atoms = 0
    for i in range(n):
        line = f.readline()
        try:
            if line[31] != 'H':
                number_of_atoms += 1
        except Exception as error:
            print("Some error occurred", line, error)

    return number_of_atoms


def dihedral(p):
    """https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python"""
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array([v - (v.dot(b[1]) / b[1].dot(b[1])) * b[1] for v in [b[0], b[2]]])
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1, 1)
    b1 = b[1] / np.linalg.norm(b[1])
    x = np.dot(v[0], v[1])
    m = np.cross(v[0], b1)
    y = np.dot(m, v[1])
    return np.degrees(np.arctan2(y, x))


def collect_atoms(molecule, ring_element):
    ring_atoms = []
    for atom_i in ring_element:
        atom = molecule.get_atom(atom_i)
        ring_atoms.append(atom)
    atom = molecule.get_atom(6)
    ring_atoms.append(atom)
    atom = molecule.get_atom(9)
    ring_atoms.append(atom)
    return ring_atoms


def calculate_bonds(ring_atoms, at_number):
    try:
        bond = ring_atoms[at_number].distance_to(ring_atoms[at_number + 1])
    except:
        bond = ring_atoms[at_number].distance_to(ring_atoms[0])
    return bond


def calculate_reference_bonds(ligand, ring_element):
    ideal_atoms = collect_atoms(ligand, ring_element)
    ideal_bonds = []
    for i in range(len(ideal_atoms)):
        ideal_bonds.append(calculate_bonds(ideal_atoms, i))
    return ideal_bonds


def average(lst):
    return sum(lst) / len(lst)


def conf_check(wings):
    if wings > 0:
        conf = "boat"
    elif wings < 0:
        conf = "chair"
    else:
        conf = "Error! Flat"
    return conf


def bonds_statistics(number_of_atoms, ligand_name, conformation_name, angle_twist, sigma_dev):
    ligand_bonds = []
    ligand_bonds_err = []
    ideal_bonds = ""
    twist_err = ""
    if conformation_name == "chair":
        ideal_bonds, twist_err = ideal_ligand_geometry(angle_twist, ligand_name, '.cif')
    elif conformation_name == "boat":
        ideal_bonds1, twist_err1 = ideal_ligand_geometry(angle_twist, ligand_name, 'cif')
        ideal_bonds2, twist_err2 = ideal_ligand_geometry(angle_twist, ligand_name, '.cif')
        if min(abs(twist_err2), abs(twist_err1)) == twist_err1:
            ideal_bonds = ideal_bonds1
            twist_err = twist_err1
        else:
            ideal_bonds = ideal_bonds2
            twist_err = twist_err2
    if np.degrees(abs(twist_err)) > sigma_dev:
        conformation_name = "twist_" + conformation_name
    for i in range(len(number_of_atoms)):
        ligand_bonds.append(calculate_bonds(number_of_atoms, i))
    for b in range(len(ligand_bonds)):
        try:
            ligand_bonds_err.append(ligand_bonds[b] - ideal_bonds[b])
        except Exception as error:
            ligand_bonds_err.append(ligand_bonds[b] - ideal_bonds[-1])
    return ligand_bonds, ligand_bonds_err, twist_err, conformation_name


def ideal_ligand_geometry(angle_twist, ligand_name, filename):
    ideal_structure = MonomerStructure.from_cif(ligand_name + filename)
    ideal_bonds = calculate_reference_bonds(ideal_structure, find_rings(ideal_structure)[0])
    ideal_ring = find_rings(ideal_structure)[0]
    ideal_atoms = collect_atoms(ideal_structure, ideal_ring)

    r = Residue(0, "ALA")
    c = Chain("A")
    r.owner(c)
    s = Structure("1xxx")
    c.owner(s)
    for a in ideal_atoms:
        a.owner(r)
    f = SaturatedRing6Geometry(ideal_atoms[0], ideal_atoms[1], ideal_atoms[2], ideal_atoms[3], ideal_atoms[4],
                               ideal_atoms[5])
    ideal_twist = f.twist_angle()
    twist_err = angle_twist - ideal_twist
    return ideal_bonds, twist_err


def atoms_bfactor(bioshell_residue):
    atoms_bf = []
    for atom in range(bioshell_residue.count_atoms()):
        atom_bf = bioshell_residue.get_atom(atom).b_factor()  # residue_type().n_atoms):
        atoms_bf.append("{:3.2f}".format(float(atom_bf)))
    return float(min(atoms_bf)), float(max(atoms_bf)), float(round(avg(atoms_bf), 2))


def avg(lst):
    return reduce(lambda a, b: float(a) + float(b), lst) / len(lst)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process geometry of ligands with six-membered ring from PDB structure.')
    parser.add_argument('pdb_model', nargs='?', default='./7ol5.pdb',
                        help='A PDB file that contains structure with the ligand of interest')
    parser.add_argument('ligand_id', nargs='?', default='EPE', type=str,
                        help='The three-letter PDB code of the ligand of interest from the given PDB file')
    parser.add_argument('-d', '--distance', nargs='?', type=int, default=5,
                        help='Radius around a ligand used to visualize the structure of the ligand and its surroundings')
    parser.add_argument('-o', '--out', nargs='?', type=argparse.FileType('w'),
                        help='The three-letter PDB code of the ligand of interest from the given PDB file')

    args = parser.parse_args()
    np.set_printoptions(precision=2)
    code = args.ligand_id
    pdb_file = args.pdb_model
    distance = 5
    ligand_in_file, fname = extract_ligand(pdb_file, code, distance)
    chain_name = ligand_in_file.back().chain_id
    res_id = ligand_in_file.id()
    chain_sel = ChainSelector(chain_name)
    ligand_sel = SelectResidueByName(code)
    filter_ligand = SelectChainResidues(chain_sel, ligand_sel)
    mol = PdbMolecule.from_pdb(fname, filter_ligand)
    ideal = download_ideal_cif('EPE')
    n_atoms = ideal.n_heavy_atoms
    sigma = 10.0  # based on https://doi.org/10.1016/j.str.2021.02.004 a 10 deg is used as deviation
    #fname = sys.argv[1]
    chain_name = fname.split("-")[2]
    code = fname.split("-")[0]
    res_id = fname.split("-")[1]
    pdb_id = fname.split("-")[3][0:4]
    #n_atoms = load_atoms_counts(PATH_TO_IDEAL_SDF + code + "_ideal.sdf")
    chain_sel = ChainSelector(chain_name)
    ligand_sel = SelectResidueByName(code)
    filter_ligand = SelectChainResidues(chain_sel, ligand_sel)
    mol = PdbMolecule.from_pdb(fname, filter_ligand)
    if mol.count_atoms() != n_atoms:
        print("Too few atoms in %s! Expected: %d Found: %d" % (fname, n_atoms, mol.count_atoms()), file=sys.stderr)
        with open(missing_atoms_file, "a") as missing_file:
            print("{:4s} {:4s} {:2s} {:3s} {:3s} {:3d} {:3d}".format(
                fname, pdb_id, chain_name, res_id, code, n_atoms, mol.count_atoms()), file=missing_file)
        sys.exit(0)

    rings = find_rings(mol)
    if len(rings) == 0:
        print("No rings found in", fname, file=sys.stderr)

    bf_min, bf_max, bf_avg = atoms_bfactor(mol)

    for ring in rings:
        if len(ring) != 6:
            print("Ring length differs from 6", fname, file=sys.stderr)
            sys.exit(0)
        atoms = collect_atoms(mol, ring)

        points = [np.array([atoms[i].x, atoms[i].y, atoms[i].z]) for i in range(8)]

        t1 = dihedral(np.array(points[0:4]))
        t2 = dihedral(np.array(points[1:5]))
        t3 = dihedral(np.array(points[2:6]))
        points = [np.array([atoms[i].x, atoms[i].y, atoms[i].z]) for i in [1, 2, 3, 6]]
        t4 = dihedral(np.array(points[0:4]))
        points = [np.array([atoms[i].x, atoms[i].y, atoms[i].z]) for i in [2, 1, 0, 7]]
        t5 = dihedral(np.array(points[0:4]))
        g = SaturatedRing6Geometry(atoms[0], atoms[1], atoms[2], atoms[3], atoms[4], atoms[5])

        w1_w2 = g.first_wing_angle() * g.second_wing_angle()
        twist_angle = g.twist_angle()
        conformation = conf_check(w1_w2)
        bonds, bonds_err, twist_angle_err, conformation = bonds_statistics(atoms, code, conformation, twist_angle,
                                                                           sigma)
        print("pdb_id", "chain", "res_no", "ligand", "conformation", "wing_atom1", "wing_atom2",
              "twist_angle", "twist_err", "t1", "t2", "t3", "avg_bond_err", "bf_min", "bf_max", "bf_avg")

        print(pdb_id, chain_name, res_id, code, conformation, g.first_wing().atom_name(),
              g.second_wing().atom_name(), round(np.degrees(twist_angle),3), round(np.degrees(twist_angle_err),3),
              round(t1,2), round(t2,2), round(t3,2), round(average(bonds_err)), bf_min, bf_max, bf_avg)
