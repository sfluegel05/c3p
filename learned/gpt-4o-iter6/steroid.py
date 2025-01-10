"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    Steroids have a modified cyclopenta[a]phenanthrene carbon skeleton.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    # Improved SMARTS for recognizing the steroid core structure and flexibility 
    # Core cyclopenta[a]phenanthrene with allowances for variations
    steroid_core_smarts = "C1CC2CC3(C)C4CCC(C3)C2C1"  # More flexible representation
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES format."

    # Create common functional group patterns for further verification
    common_steroid_features = [
        Chem.MolFromSmarts("C=O"),
        Chem.MolFromSmarts("O"),
        Chem.MolFromSmarts("C[C@H](O)C(=O)"),
        Chem.MolFromSmarts("CC(C)(C)O"),
        Chem.MolFromSmarts("CC(=O)O"),
        Chem.MolFromSmarts("C#CC[C@H]"),
    ]
    
    # Check if the molecule matches the steroid core
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Does not contain the flexible steroid core structure."

    # Evaluate the presence of typical steroid modifications
    for feature in common_steroid_features:
        if mol.HasSubstructMatch(feature):
            return True, "Contains characteristic functional groups common in steroids."

    return True, "Matches flexible steroid core."

# Some example SMILES for test
examples = [
    "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)[C@H](F)C[C@]34C)[C@@H]1CC[C@@H]2O",  # Steroid
    "CC(=O)O[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CCCC[C@]4(C)[C@H]3CC[C@]12C"         # Steroid
]

for smiles in examples:
    result, reason = is_steroid(smiles)
    print(f"SMILES: {smiles} -> Is Steroid: {result}, Reason: {reason}")