"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    
    # Define a general inositol pattern capturing potential stereochemistry variations
    inositol_pattern = Chem.MolFromSmarts("OC[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No 1D-myo-inositol moiety found"

    # Define a pattern for detecting phosphate ester linkage to inositol
    phosphate_linkage_pattern = Chem.MolFromSmarts("O[P](=O)(O)OC[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(phosphate_linkage_pattern):
        return False, "No phosphate ester linkage to inositol at position 1"
    
    # Define pattern for glycerol backbone in a phosphatidyl linkage
    glycerol_pattern = Chem.MolFromSmarts("O[C@@H]COC(=O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    return True, "Contains 1-phosphatidyl-1D-myo-inositol moiety"