"""
Classifies: CHEBI:74716 withanolide
"""
from rdkit import Chem

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is defined as a C28 steroid lactone with modified side chains forming a lactone ring.

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
    
    # Check for steroid core (4 fused rings pattern, allowing modification)
    steroid_patterns = [
        # Pattern for generic steroid core structure
        Chem.MolFromSmarts("C1CC2CCC3C4CCC(C4)C3C2C1"),
        Chem.MolFromSmarts("C1CCCC2C1C3CCC4C3C2CC4")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in steroid_patterns):
        return False, "No steroid core structure found"

    # Check for presence of lactone ring
    lactone_pattern = Chem.MolFromSmarts("O=C1OC=CC1")
    if not mol.HasSubstructMatch(lactone_pattern):
        return False, "No lactone ring found"

    # Check for C28 steroid lactone
    carbon_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if carbon_count != 28:
        return False, f"Molecule has {carbon_count} carbons, expected 28"

    # Additional checks for common withanolide functional groups
    hydroxyl_pattern = Chem.MolFromSmarts("C(C)([OH])[C]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No common withanolide hydroxyl groups found"

    return True, "The molecule contains key features of a withanolide"

# Example usage
example_smiles = "CC1=C(C)C(=O)O[C@H](C1)[C@](C)(O)[C@@]12O[C@@H]1C[C@H]1[C@@H]3C[C@H]4O[C@]44[C@@H](O)C=CC(=O)[C@]4(C)[C@H]3CC[C@]21C"
result = is_withanolide(example_smiles)
print(result)