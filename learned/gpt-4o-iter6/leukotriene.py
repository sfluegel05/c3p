"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    Leukotrienes have a characteristic polyunsaturated C20 fatty acid backbone
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

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, f"Contains {c_count} carbon atoms, expected at least 20"

    # Identify double bonds and conjugation
    db_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if db_count < 4:
        return False, f"Contains {db_count} double bonds, expected at least 4"

    # Check for conjugated double bonds
    # Pattern for three conjugated double bonds (with potential variations)
    conjugated_pattern = Chem.MolFromSmarts('C=CC=CC=C')
    if not mol.HasSubstructMatch(conjugated_pattern):
        # Consider alternative patterns such as presence of cycles in conjugation
        conjugated_cycle_pattern = Chem.MolFromSmarts('C=C(-C=C)-C=C')
        if not mol.HasSubstructMatch(conjugated_cycle_pattern):
            return False, "Does not contain at least three conjugated double bonds"

    # Check for additional known functional groups or patterns
    # Example: Terminal carboxylic acid group often present in leukotrienes
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)[O-]') 
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "Lacks typical terminal carboxylic group found in leukotrienes"

    return True, "Contains characteristic C20 backbone with conjugated double bond pattern and known functional groups"