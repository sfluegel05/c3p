"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define SMARTS patterns for the structural features of 1-phosphatidyl-1D-myo-inositol
    # Glycerol backbone with ester linkages
    glycerol_ester_pattern = Chem.MolFromSmarts("[C@H](OC(=O))COC(=O)")
    if not mol.HasSubstructMatch(glycerol_ester_pattern):
        return False, "No glycerol backbone with ester linkages found"

    # Phosphate group pattern: phosphate group attached to the glycerol
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group attached to the glycerol backbone"

    # Inositol ring pattern: cyclohexane ring with hydroxyls
    inositol_pattern = Chem.MolFromSmarts("C1(C(O)C(O)C(O)C(O)C1O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol cyclohexane ring with hydroxyl groups found"

    return True, "Matches 1-phosphatidyl-1D-myo-inositol structural features"

# Example usage:
smiles_examples = [
    "CCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)OC(=O)CCCCCCCCCCC",
    "CCCCCCCCC\C=C/C\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)OC(=O)CCCCCCC"
]

for example in smiles_examples:
    result, reason = is_1_phosphatidyl_1D_myo_inositol(example)
    print(f"SMILES: {example} -> is 1-phosphatidyl-1D-myo-inositol: {result}, Reason: {reason}")