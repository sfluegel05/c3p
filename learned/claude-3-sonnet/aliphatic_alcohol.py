"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: CHEBI:15756 aliphatic alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is an alcohol derived from an aliphatic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for aliphatic carbon chains (including branched and unsaturated)
    aliphatic_pattern = Chem.MolFromSmarts("[CX4H3,CX4H2,CX3]~[CX4H2,CX3]~[CX4H2,CX3]")
    aliphatic_matches = mol.GetSubstructMatches(aliphatic_pattern)
    
    if not aliphatic_matches:
        return False, "No aliphatic carbon chain found"

    # Look for alcohol functional group (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OX1H]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    
    if not alcohol_matches:
        return False, "No alcohol functional group found"

    # Check if the alcohol group is attached to the aliphatic chain
    for alcohol_match in alcohol_matches:
        atom = mol.GetAtomWithIdx(alcohol_match[0])
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic() == False:
                break
        else:
            continue  # Alcohol group not attached to aliphatic chain
        return True, "Contains an aliphatic carbon chain with an alcohol functional group"

    return False, "Alcohol group not attached to aliphatic carbon chain"