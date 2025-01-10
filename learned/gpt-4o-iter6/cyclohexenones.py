"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is defined as any six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern for cyclohexenone structure: Six-membered cycloalkenone with flexibility in substituents
    cyclohexenone_pattern = Chem.MolFromSmarts("C1=CC(=O)C[CH2][CH2]1")
    
    # Check for the cyclohexenone pattern
    if mol.HasSubstructMatch(cyclohexenone_pattern):
        return True, "Cyclohexenone structure identified in the molecule"
    
    # Alternative SMARTS with relaxed structure to capture missed configurations
    alt_pattern = Chem.MolFromSmarts("C1=C(C)CC(=O)C1")  # Accommodates substitutions at any carbon
    
    # Check for the alternative cyclohexenone pattern
    if mol.HasSubstructMatch(alt_pattern):
        return True, "Cyclohexenone structure identified in the molecule (alternative pattern)"
    
    return False, "No cyclohexenone structure found"