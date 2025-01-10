"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdqueries

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is any molecule that contains two amino-acid residues connected by peptide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dipeptide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define alpha-amino acid residue pattern
    amino_acid_pattern = Chem.MolFromSmarts("""
        [NX3,NX4+][CX4H]([*])[CX3](=O)[O-,OH]
        |
        [NX3,NX4+][CX4H]([*])[CX3](=O)[O-,OH]
    """)

    # Define peptide bond pattern (excluding side-chain amides)
    peptide_bond_pattern = Chem.MolFromSmarts("[$([CX3](=O)[NX3H])]")

    # Find all peptide bonds
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    n_peptide_bonds = len(peptide_bonds)

    # Check for at least one peptide bond
    if n_peptide_bonds == 0:
        return False, "No peptide bonds found"

    # Attempt to identify amino acid residues
    # We will label the atoms to keep track of residues
    atom_labels = {}
    for atom in mol.GetAtoms():
        atom.SetProp('residue', '-1')  # initialize with -1

    residue_idx = 0
    for match in mol.GetSubstructMatches(amino_acid_pattern):
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            atom.SetProp('residue', str(residue_idx))
        residue_idx += 1

    n_residues = residue_idx

    if n_residues != 2:
        return False, f"Found {n_residues} amino acid residues, expected 2 for a dipeptide"

    # Check connectivity between residues via peptide bonds
    # Get the residue labels for atoms involved in peptide bonds
    residues_connected = set()
    for match in peptide_bonds:
        c_idx = match[0]
        n_idx = match[1]
        c_atom = mol.GetAtomWithIdx(c_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)
        res_c = c_atom.GetProp('residue')
        res_n = n_atom.GetProp('residue')
        if res_c != res_n:
            residues_connected.update([res_c, res_n])

    if len(residues_connected) != 2:
        return False, "Amino acid residues are not properly connected via peptide bonds"

    return True, "Molecule is a dipeptide composed of two amino acid residues connected by peptide bonds"