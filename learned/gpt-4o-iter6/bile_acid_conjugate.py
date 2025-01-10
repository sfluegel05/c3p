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

    # Generalized steroid backbone pattern to cover bile acids' variety
    steroid_backbone_pattern = Chem.MolFromSmarts("C1(C2CCC3C(C(O)CCC3)C2)CCCCC1")  # More inclusive pattern
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No bile acid core structure found"

    # Expanded conjugate patterns with flexibility
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)[O-]")  # Carboxylate form
    taurine_pattern = Chem.MolFromSmarts("S(=O)(=O)NCC(=O)")  # More inclusive sulfonamide group
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O-]")  # Charged sulfate group
    glucuronic_acid_pattern = Chem.MolFromSmarts("OC[C@@H](O)[C@@H](O)C(=O)[O-]")  # Flexible ester linkage

    # Check for presence of any known conjugation pattern
    conjugates_found = False
    if mol.HasSubstructMatch(glycine_pattern):
        conjugates_found = True
    elif mol.HasSubstructMatch(taurine_pattern):
        conjugates_found = True
    elif mol.HasSubstructMatch(sulfate_pattern):
        conjugates_found = True
    elif mol.HasSubstructMatch(glucuronic_acid_pattern):
        conjugates_found = True

    if not conjugates_found:
        return False, "No known bile acid conjugate pattern found"

    return True, "Matches bile acid core with conjugate pattern"

# Test with example SMILES
example_smiles = "C1[C@@]2(CCC(=O)NCC(=O)O)C(C)CCC3C4(C)CCC(=O)CCC4CCC23"
result, reason = is_bile_acid_conjugate(example_smiles)
print(result, reason)