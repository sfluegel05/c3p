"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:27270 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is an organic amino compound with an amino group and a hydroxy group
    attached to the same carbon atom, which may or may not have additional substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for atoms with both amino (-NH2/-NH-) and hydroxy (-OH) groups attached to the same carbon
    hemiaminal_pattern = Chem.MolFromSmarts("[NH2,NH1][CX4;H1]([OH1])")
    hemiaminal_atoms = mol.GetSubstructMatches(hemiaminal_pattern)
    
    if hemiaminal_atoms:
        return True, "Contains amino and hydroxy groups attached to the same carbon atom"
    else:
        return False, "Does not contain amino and hydroxy groups attached to the same carbon atom"