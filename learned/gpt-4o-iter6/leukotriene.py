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

    # Count carbon atoms to ensure C20 backbone
    c_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if c_count != 20:
        return False, f"Contains {c_count} carbon atoms, expected 20 for leukotrienes"

    # Identify double bonds, ensuring four are present
    db_count = sum(bond.GetBondType() == Chem.BondType.DOUBLE for bond in mol.GetBonds())
    if db_count < 4:
        return False, f"Contains {db_count} double bonds, expected at least 4"

    # Improved check for three conjugated double bonds
    conjugated_pattern = Chem.MolFromSmarts('[#6]=[#6]-[#6]=[#6]-[#6]=[#6]')
    conjugated_patterns_variants = [
        '[#6]=[#6]-[#6]=[#6]-[#6]=[#6]', 
        '[#6]=[#6]~[#6]=[#6]-[#6]=[#6]', 
        '[#6]C=C=C=CC=[#6]'
    ]
    if not any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in conjugated_patterns_variants):
        return False, "Fails to have three conjugated double bonds"

    # Check for terminal carboxylic acid group or variants
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)O')  # Match carboxylic acids including charged forms
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Lacks terminal carboxylic acid group"

    return True, "Contains characteristic C20 backbone with conjugated double bond pattern and carboxylic group"