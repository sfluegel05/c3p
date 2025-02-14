"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:18154 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define amino group SMARTS pattern (primary and secondary amines)
    # Match nitrogen atoms that are sp3 hybridized, not aromatic, not double-bonded, and with at least one hydrogen
    amino_smarts = "[#7X3&H1,#7X3&H2]"  # Nitrogen with one or two hydrogens attached
    amino_pattern = Chem.MolFromSmarts(amino_smarts)
    if amino_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Find all amino groups in the molecule
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    num_amino_groups = len(amino_matches)

    if num_amino_groups >= 2:
        return True, f"Molecule contains {num_amino_groups} amino groups"
    else:
        return False, f"Molecule contains {num_amino_groups} amino group(s), need at least 2"