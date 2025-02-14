"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for zwitterionic alpha-amino acids
    # Look for alpha-amino group and deprotonated carboxylate as core structure
    core_zwitterion_pattern = Chem.MolFromSmarts('[C@@H]([NH3+])C(=O)[O-]')
    
    # Complex side chain matching with flexibility for backbone
    full_zwitterion_patterns = [
        Chem.MolFromSmarts('[$([NH3+]),$([NH2+]),$([NH+]=*)][C@H](C[*])C(=O)[O-]'),
        Chem.MolFromSmarts('[C@@H]([NH3+])(C[*])C(=O)[O-]'),
    ]

    # Check for presence of zwitterionic core
    if not mol.HasSubstructMatch(core_zwitterion_pattern):
        return False, "No matching alpha-amino acid zwitterion core structure found"
    
    # Check each potential pattern to ensure there's commonly encountered complexity
    for pattern in full_zwitterion_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains alpha-amino acid structure with zwitterionic charges"

    return False, "Missing patterns indicating complex zwitterionic alpha-amino acid configuration"

# Example usage:
# smiles = "CC(=O)CC([NH3+])C([O-])=O"  # Example SMILES from the list
# result, reason = is_alpha_amino_acid_zwitterion(smiles)
# print(f"Result: {result}, Reason: {reason}")