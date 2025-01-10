"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is defined as a furan ring fused with a coumarin structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a furanocoumarin substructure pattern 
    # This pattern attempts to identify key structural elements potentially shared by examples
    # Fused 5-membered furan ring to a 6-membered lactone ring
    furanocoumarin_pattern = Chem.MolFromSmarts('C1OC=2C=CC(=O)OC2=C1')

    # Check for the presence of the furanocoumarin fundamental fused ring system
    matches = mol.GetSubstructMatches(furanocoumarin_pattern)
    if not matches:
        return False, "No furanocoumarin core structure found"

    return True, "Contains furanocoumarin core fused ring structure"

# This function checks specifically for the fused structure indicative of a furanocoumarin. 
# Depending on SMILES complexity encountered in specific instances, further refinement of SMARTS patterns and chemical logic may be required.