"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has a hydroxy group attached to a carbon which is itself bonded to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for an aromatic carbon directly bonded to a hydroxyl group
    primary_alcohol_pattern = Chem.MolFromSmarts("[OHX2][CX4]")  
    aromatic_carbon_pattern = Chem.MolFromSmarts("[#6r][OH]")  # aromatic carbon bonded to hydroxyl group

    # Check if the molecule matches primary alcohol pattern
    if not mol.HasSubstructMatch(primary_alcohol_pattern):
        return False, "No primary alcohol (C-OH) group detected"

    # Check if the molecule contains an aromatic ring
    if not Chem.GetSSSR(Chem.MolFromSmarts("a") & mol):
        return False, "No aromatic ring detected"
    
    # Check for primary alcohol group being attached to an aromatic carbon
    if not mol.HasSubstructMatch(aromatic_carbon_pattern):
        return False, "Alcohol group not attached to aromatic carbon"

    return True, "Contains hydroxy group attached to a carbon bonded to an aromatic ring"

# The function can now be used to classify SMILES strings for aromatic primary alcohols.