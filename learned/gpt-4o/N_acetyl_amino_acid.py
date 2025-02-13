"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid has an acetyl group linked to the nitrogen of an amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an N-acetyl group specifically bound to a nitrogen atom
    n_acetyl_pattern = Chem.MolFromSmarts("CN(C=O)")

    # Check if the N-acetyl pattern is present
    if not mol.HasSubstructMatch(n_acetyl_pattern):
        return False, "No N-acetyl group directly on nitrogen found"

    # Look for common amino acid backbone or analog structures
    amino_acid_patterns = [
        Chem.MolFromSmarts("N[C@@H](C(=O)O)"),   # L-amino acid pattern
        Chem.MolFromSmarts("N[C@H](C(=O)O)"),    # D-amino acid pattern
        Chem.MolFromSmarts("NCC(=O)O"),          # Generic amino acid pattern
        Chem.MolFromSmarts("N1CCCC1C(=O)O"),     # Cyclic (proline-like) pattern
        Chem.MolFromSmarts("NC(=O)[C@H](C(=O)O)"),  # Variants with additional attached groups
    ]
    
    has_amino_acid_structure = any(mol.HasSubstructMatch(pattern) for pattern in amino_acid_patterns)
    if not has_amino_acid_structure:
        return False, "No identifiable amino acid structure"

    # Allow for moderate molecular size to accommodate known examples
    if mol.GetNumAtoms() > 35:  # Updated limit
        return False, "Molecule too large, likely not a simple N-acetyl-amino acid"

    return True, "Contains N-acetyl group attached to an amino acid structure"