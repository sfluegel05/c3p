"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester has a decanoyl group (10 carbons) attached via an ester bond.

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

    # Define the SMARTS pattern for the ester group
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")

    # Find all ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester bond found"
    
    # Check for the decanoyl group (10 carbon chain attached to C=O) for each ester
    decanoyl_pattern = Chem.MolFromSmarts("[CX3](=O)[CX4]([CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4])")
    
    for match in ester_matches:
          
        # Now check for the decanoyl group as a substructure attached to the ester group
        if mol.HasSubstructMatch(decanoyl_pattern):

             return True, "Contains a decanoyl group connected via an ester bond"
    
    
    return False, "No decanoate ester core found"