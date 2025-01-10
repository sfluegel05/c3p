"""
Classifies: CHEBI:86315 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide where at least one of the groups attached to sulfur is a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the methyl sulfide pattern (sulfur bonded to a methyl group)
    methyl_sulfide_pattern = Chem.MolFromSmarts("CS")
    
    # Check for the methyl sulfide pattern
    if mol.HasSubstructMatch(methyl_sulfide_pattern):
        return True, "Contains a sulfur atom bonded to a methyl group"
    
    return False, "Does not contain a sulfur atom bonded to a methyl group"

# Example usage:
# smiles = "CCSC"  # Example SMILES for ethyl methyl sulfide
# result, reason = is_methyl_sulfide(smiles)
# print(result, reason)