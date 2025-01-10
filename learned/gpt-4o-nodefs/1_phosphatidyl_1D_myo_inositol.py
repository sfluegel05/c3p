"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a given SMILES string corresponds to a chemical that is a 1-phosphatidyl-1D-myo-inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glycerol backbone with two ester linkages and a phospho group
    glycerol_phospho_pattern = Chem.MolFromSmarts("O[C@@H]1COC2(=O)COP(=O)(O)O2")
    if not mol.HasSubstructMatch(glycerol_phospho_pattern):
        return False, "No glycerol backbone with phosphatidyl linkages found"

    # Define the myo-inositol ring structure with specific stereochemistry
    myo_inositol_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol ring found with correct hydroxyl configurations"

    return True, "Successfully matches 1-phosphatidyl-1D-myo-inositol structural features"

# Example usage:
smiles_examples = [
    "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)OC(=O)CCCCCCCCCCC",
    "CCCCCCCCC\C=C/C\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)OC(=O)CCCCCCC"
]

for example in smiles_examples:
    result, reason = is_1_phosphatidyl_1D_myo_inositol(example)
    print(f"SMILES: {example} -> is 1-phosphatidyl-1D-myo-inositol: {result}, Reason: {reason}")