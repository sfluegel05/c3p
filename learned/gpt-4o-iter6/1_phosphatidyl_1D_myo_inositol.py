"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    The molecule should have a 1D-myo-inositol with a phosphatidyl group at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-phosphatidyl-1D-myo-inositol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Correctly map 1D-myo-inositol based on stereochemistry
    inositol_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No 1D-myo-inositol moiety found"

    # Phosphate group attached at position 1 of inositol
    phosphate_linkage_pattern = Chem.MolFromSmarts("O[P](=O)(O)O[C@H]1")
    if not mol.HasSubstructMatch(phosphate_linkage_pattern):
        return False, "No phosphate ester linkage at position 1 of inositol"

    # Detect glycerol backbone linked to the phosphate (phosphatidyl group)
    glycerol_pattern = Chem.MolFromSmarts("OC[C@H]COC(=O)C")  # Minimal glycerol scaffold in phosphatidyl linkage
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone for phosphatidyl linkage"

    return True, "Contains 1-phosphatidyl-1D-myo-inositol moiety"