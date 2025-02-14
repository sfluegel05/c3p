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
    # Look for protonated amine on alpha carbon and a deprotonated carboxylate
    alpha_amino_acid_pattern = Chem.MolFromSmarts('[$([NH3+]),$([NH2+]),$([NH+]=*)]CC(=O)[O-]')

    if not mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return False, "No matching alpha-amino acid zwitterion pattern found"

    return True, "Contains alpha-amino acid structure with zwitterionic charges"

# Example usage:
# smiles = "CC(=O)CC([NH3+])C([O-])=O"  # Example SMILES from the list
# result, reason = is_alpha_amino_acid_zwitterion(smiles)
# print(f"Result: {result}, Reason: {reason}")