"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    Withanolides are C28 steroid lactones with modified side chains forming lactone rings and substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for the basic steroidal backbone (four fused rings)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CC3C(C2)CCCC3")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone detected"
    
    # Check for the lactone ring (cyclic ester)
    lactone_pattern = Chem.MolFromSmarts("C1OC(=O)C=CC1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone group found"
    
    # Count the carbon atoms to check for a C28 steroid
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 28:
        return False, f"Not enough carbon atoms: found {c_count}, expected at least 28"

    # Check for hydroxyl and ketone functional groups
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    if not hydroxyl_matches and not ketone_matches:
        return False, "No hydroxyl or ketone functionalities detected"

    return True, "Contains features consistent with withanolides (steroid backbone, lactone group, C28)"