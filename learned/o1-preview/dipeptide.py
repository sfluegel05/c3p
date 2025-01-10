"""
Classifies: CHEBI:46761 dipeptide
"""
"""
Classifies: dipeptide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dipeptide(smiles: str):
    """
    Determines if a molecule is a dipeptide based on its SMILES string.
    A dipeptide is any molecule that contains two amino-acid residues connected by a peptide linkage.

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

    # Define peptide bond pattern (amide bond between C=O and N-H)
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3][CX4]")
    if peptide_bond_pattern is None:
        return False, "Failed to create peptide bond SMARTS pattern"

    # Find all peptide bonds
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    n_peptide_bonds = len(peptide_bonds)

    # Check for exactly one peptide bond
    if n_peptide_bonds != 1:
        return False, f"Found {n_peptide_bonds} peptide bonds, expected exactly 1 for a dipeptide"

    # Get the bond indices of the peptide bond(s)
    peptide_bond_indices = []
    for match in mol.GetSubstructMatches(peptide_bond_pattern, uniquify=True):
        # Get the indices of the carbonyl carbon and the amide nitrogen
        c_idx = match[0]
        n_idx = match[1]
        # Find the bond between these atoms
        bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
        if bond is not None:
            peptide_bond_indices.append(bond.GetIdx())

    if len(peptide_bond_indices) != 1:
        return False, "Could not uniquely identify the peptide bond"

    # Fragment the molecule at the peptide bond(s)
    fragmented_mol = Chem.FragmentOnBonds(mol, peptide_bond_indices, addDummies=True)
    fragments = Chem.GetMolFrags(fragmented_mol, asMols=True)

    n_fragments = len(fragments)

    # Check for exactly two fragments (amino acid residues)
    if n_fragments != 2:
        return False, f"Found {n_fragments} fragments after breaking peptide bond, expected 2"

    # Optional: Validate that each fragment resembles an amino acid residue
    amino_acid_residue_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4][CX3](=O)")
    if amino_acid_residue_pattern is None:
        return False, "Failed to create amino acid residue SMARTS pattern"

    for idx, frag in enumerate(fragments):
        if not frag.HasSubstructMatch(amino_acid_residue_pattern):
            return False, f"Fragment {idx+1} does not resemble an amino acid residue"

    return True, "Molecule is a dipeptide composed of two amino acid residues connected by a peptide bond"