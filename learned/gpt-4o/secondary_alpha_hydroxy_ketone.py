"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone (acyloin) has a carbonyl group and hydroxy group linked by
    a secondary carbon bearing one hydrogen and one organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Detect a secondary carbon with both a hydroxyl and a carbonyl group
    # [CH] - secondary carbon, [CX3](=O) - carbonyl, [OX2H] - hydroxy connected at the alpha position
    secondary_alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts("[CH1]([OX2H])[CX3](=O)")
    
    # Check for the presence of the substructure pattern
    if mol.HasSubstructMatch(secondary_alpha_hydroxy_ketone_pattern):
        return True, "Molecule contains the characteristic structure of a secondary alpha-hydroxy ketone"
    
    return False, "No characteristic secondary alpha-hydroxy ketone structure found"

# Example test with a known secondary alpha-hydroxy ketone
test_smiles = "OC(C(O)=O)C(=O)CCC(O)=O"  # Example of secondary alpha-hydroxy ketone
print(is_secondary_alpha_hydroxy_ketone(test_smiles))