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

    # More flexible bile acid steroid backbone pattern
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CC2CCC3C(C2C1)C4CCC5C(C=CCC5C4)C3C")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No bile acid core structure found"

    # Improved conjugation patterns
    glycine_pattern = Chem.MolFromSmarts("NC[C@@H](C(=O)O)")
    taurine_pattern = Chem.MolFromSmarts("S(=O)(=O)CCN")
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[O-]")  # Includes common charged state
    generic_amino_acid_pattern = Chem.MolFromSmarts("NC(C(=O)O)C")  # Covers more generic amino acids

    # Check for presence of any known conjugation pattern
    conjugates_found = False
    if mol.HasSubstructMatch(glycine_pattern):
        conjugates_found = True
    elif mol.HasSubstructMatch(taurine_pattern):
        conjugates_found = True
    elif mol.HasSubstructMatch(sulfate_pattern):
        conjugates_found = True
    elif mol.HasSubstructMatch(generic_amino_acid_pattern):
        conjugates_found = True

    if not conjugates_found:
        return False, "No known bile acid conjugate pattern found"

    return True, "Matches bile acid core with conjugate pattern"

# Example SMILES for testing
example_smiles = "C[C@H](CCC(=O)NCC(O)=O)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C[C@@H]4C[C@H](O)CC[C@]4(C)[C@H]3C[C@H](O)[C@]12C"  # glycocholic acid
result, reason = is_bile_acid_conjugate(example_smiles)
print(result, reason)