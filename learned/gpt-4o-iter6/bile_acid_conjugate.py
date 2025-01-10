"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate involves a bile acid structure attached to hydrophilic/charged 
    groups such as glycine, taurine, sulfate, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generic bile acid steroid backbone pattern (allows for flexibility)
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CCC3C2CC4C(C3)CC(C4)C")  # Core steroid structure
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No bile acid core structure found"

    # Enhanced conjugation patterns
    conjugation_patterns = [
        Chem.MolFromSmarts("NC(C)=O"),  # Simple glycine pattern
        Chem.MolFromSmarts("S(=O)(=O)CCN"),  # Taurine pattern
        Chem.MolFromSmarts("OS(=O)(=O)O"),  # Sulfate pattern
        Chem.MolFromSmarts("C(O)C(=O)O"),  # Generic carboxylic acid pattern for glucuronates
        Chem.MolFromSmarts("C1OC(CO)C(O)C1O"),  # Glucose-like pattern
        Chem.MolFromSmarts("O[C@@H]1[C@H]([C@H](O)[C@H](O)[C@H](O)[C@H]1O)C"),  # Specific sugar motifs
        Chem.MolFromSmarts("NC(C(=O)O)C"),  # Generic amino acid pattern
    ]

    # Check for any known conjugation pattern
    conjugates_found = False
    for pattern in conjugation_patterns:
        if mol.HasSubstructMatch(pattern):
            conjugates_found = True
            break

    if not conjugates_found:
        return False, "No known bile acid conjugate pattern found"

    return True, "Matches bile acid core with conjugate pattern"

# Example SMILES for testing
example_smiles = "C[C@H](CCC(=O)NCC(O)=O)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C"  # glycocholic acid
result, reason = is_bile_acid_conjugate(example_smiles)
print(result, reason)