"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.

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
    
    # Flexible pattern for steroid core (a four-ring structure, considering stereochemistry)
    steroid_patterns = [
        Chem.MolFromSmarts("[#6]1-[#6]-[#6]-[#6]2-[#6]-[#6]-[#6]-3-[#6]-[#6]-[#6]-[#6]-4-[#6]-[CH2]-[CH]-1-[C](=O)2=3-4"), # More flexible steroid core
        Chem.MolFromSmarts("[#6]12[#6](#[#6])[#6]-[#6]-3=[#6]-[#6]-[#6]-[#6]-4-[#6]-[CH2]-[C]1-[#6](=O)2-[C]=C3=C4")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return False, "No valid steroid core structure found"
    
    # Check for presence of a lactone ring
    lactone_pattern = Chem.MolFromSmarts("[OH0,C](=O)[C&R1][O&R1;!c]")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Check for C28 steroid lactone (with some flexibility)
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if carbon_count < 26 or carbon_count > 30:
        return False, f"Molecule has {carbon_count} carbons, expected around 28"

    # Check for common functional groups that are present in withanolides
    hydroxyl_pattern = Chem.MolFromSmarts("[CH]([OH])[C&R1]")  # Hydroxyl group
    methoxy_pattern = Chem.MolFromSmarts("[CH3][OH]")  # Methoxy group
    if not (mol.HasSubstructMatch(hydroxyl_pattern) or mol.HasSubstructMatch(methoxy_pattern)):
        return False, "No hydroxyl or methoxy groups typically found in withanolides"

    return True, "The molecule contains key features of a withanolide"

# Example usage
example_smiles = "CC1=C(C)C(=O)O[C@H](C1)[C@](C)(O)[C@@]12O[C@@H]1C[C@H]1[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3CC[C@]21C"
result = is_withanolide(example_smiles)
print(result)