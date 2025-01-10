"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid contains a fatty acid chain with a carboxylic acid group and a ring structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)O)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for a ring structure
    if not mol.GetRingInfo().NumRings():
        return False, "No ring structure found"
    
    return True, "Contains a carboxylic acid group and a ring structure, classifying as a cyclic fatty acid"

# Example usage:
# result, reason = is_cyclic_fatty_acid("C1(C(C/C=C\CCCCCO)O1)CCCCCCCC(=O)O")
# print(result, reason)