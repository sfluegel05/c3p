"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a hydroxy group in the alpha- or 2-position relative to the carboxylic acid group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for both saturated and unsaturated carboxylic acid chains
    # Match carboxylic acid group (C(=O)O) and ensure the hydroxy is at the second carbon
    pattern = Chem.MolFromSmarts("[C;!R](O)[C;!R;!D4](O)[C;!R](=O)O")
    
    # Ensure no complex rings or other large functionalities exist
    complex_functional_groups = Chem.MolFromSmarts("[R]=[!N,#7,#8,#16]") # Any rings or non-S/O/N heteroatoms 
    if mol.HasSubstructMatch(complex_functional_groups):
        return False, "Molecule contains complex or cyclic functionalities typical of multi-ring structures"

    # Perform the matching
    if mol.HasSubstructMatch(pattern):
        return True, "Contains 2-hydroxy group with carboxylic acid at the 2-position"
    else:
        return False, "No 2-hydroxy group found at the alpha position"

# Example usage
print(is_2_hydroxy_fatty_acid("CCCCCCC[C@@H](O)C(O)=O"))  # (R)-2-hydroxynonanoic acid