"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    A leukotriene should have the characteristic polyunsaturated C20 fatty acid backbone
    with four double bonds, including three conjugated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Ensure there is at least the characteristic 20 carbon backbone
    if c_count < 20:
        return False, f"Too few carbon atoms: {c_count} (expected 20 or more with functional groups)"
    
    # Identify double bonds
    db_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    
    # Ensure at least four double bonds
    if db_count < 4:
        return False, f"Too few double bonds: {db_count} (expected 4 or more)"
    
    # Check for conjugated double bonds using SMARTS pattern for conjugation
    conjugated_pattern = Chem.MolFromSmarts('C=C-C=C-C=C')
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Does not contain at least three conjugated double bonds"

    return True, "Contains at least C20 backbone with characteristic conjugated double bond pattern"