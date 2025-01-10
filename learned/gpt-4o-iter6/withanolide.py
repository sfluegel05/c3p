"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    
    A withanolide is characterized as a steroid lactone with specific structural features,
    such as a C28 steroid with a modified side chain forming a lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Expanded steroid core patterns to include more diverse arrangements
    steroid_core_patterns = [
        Chem.MolFromSmarts("C1CC2=C(C3CC[C@@H]4C=C[C@H]4[C@@H]3C2)C=CC=C1"),  # 6-6-6-5 with conjugations
        Chem.MolFromSmarts("C1CCC2C3CCC4(C)C=CC4C3CC2C1C"),                 # Common steroid nucleus
    ]
    
    # Check if any steroid core pattern is present
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_core_patterns):
        return False, "No valid steroid core structure found"
    
    # Flexible lactone recognition
    lactone_patterns = [
        Chem.MolFromSmarts("O=C1O[C@@H]1"),            # 3 to 7-member lactone configurations
        Chem.MolFromSmarts("O=C1COC=C1"),              # possibility of substitutions
    ]
    
    # Check for lactone presence
    if not any(mol.HasSubstructMatch(pattern) for pattern in lactone_patterns):
        return False, "No lactone ring found"
    
    # Count total carbons
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if not (26 <= carbon_count <= 30):
        return False, f"Molecule has {carbon_count} carbons, outside expected range for withanolides"
    
    # Presence of hydroxyl groups matching common withanolide substituent
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4H1][OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Return True if all conditions hold
    return True, "The molecule contains key features of a withanolide"