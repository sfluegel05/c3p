"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: CHEBI:16670 polypeptide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide based on its SMILES string.
    A polypeptide is a peptide containing ten or more amino acid residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broad pattern to match amino acid residues
    # Looks for N-C(=O)-C pattern, which is common in amino acids
    aa_pattern = Chem.MolFromSmarts("[NX3H1,nX3H0][CX4H1]([CX3](=[OX1]))[CX4H1R0]")
    aa_matches = mol.GetSubstructMatches(aa_pattern)
    
    # Each match corresponds to one amino acid residue
    num_residues = len(aa_matches)
    
    if num_residues < 10:
        return False, f"Found {num_residues} amino acid residues, need at least 10"
    
    # Check for peptide bonds
    peptide_bond_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[NX3H1,nX3H0]")
    peptide_bonds = len(mol.GetSubstructMatches(peptide_bond_pattern))
    
    # Ensure there are at least (num_residues - 1) peptide bonds
    if peptide_bonds < num_residues - 1:
        return False, "Insufficient peptide bonds for polypeptide"
    
    # Check for typical polypeptide elements
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < num_residues or o_count < num_residues:
        return False, "Insufficient nitrogen or oxygen atoms for polypeptide"
    
    return True, f"Contains {num_residues} amino acid residues, qualifies as a polypeptide"