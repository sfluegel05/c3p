"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester results from the esterification of decanoic acid with an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for decanoate ester: Ensure exactly 10 carbons in chain with ester linkage
    # Account for ester groups as well as potential branching/double bonds
    decanoate_backbone_pattern = Chem.MolFromSmarts("[#6]1-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6](=O)-[#8]") # Flexible for potential modifications in chain
    
    # Check decanoate backbone exists
    if mol.HasSubstructMatch(decanoate_backbone_pattern):
        return True, "Contains decanoate ester structure: recognized as decanoate ester"
    
    return False, "No decanoate ester structure found"