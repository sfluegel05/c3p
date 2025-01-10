"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    This is identified by a glycerol backbone, acyl group, and phosphoethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Improved glycerol backbone pattern: ensuring chirality and glycerol marks accurately
    glycerol_pattern = Chem.MolFromSmarts("C(CO)[C@@H]O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Improved acyl pattern to match a carbonyl group bound to the remaining glycerol
    acyl_pattern = Chem.MolFromSmarts("C(=O)OC[C@H]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No acyl group found"
    
    # Enhanced pattern for phosphoethanolamine, accounting for possible ionization
    phospho_ethanolamine_pattern = Chem.MolFromSmarts("COP(O)(=O)OCC[N]")
    if not mol.HasSubstructMatch(phospho_ethanolamine_pattern):
        return False, "No phosphoethanolamine group found"

    # If all features are present, classify as 1-acyl-sn-glycero-3-phosphoethanolamine
    return True, "Contains glycerol backbone, acyl group, and phosphoethanolamine group"