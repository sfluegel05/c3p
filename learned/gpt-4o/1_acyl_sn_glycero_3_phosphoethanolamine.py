"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise.
        str: Reason for classification.
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for a glycerol backbone with the correct chiral center and phosphate-ethanolamine linkage
    glycerol_pattern = Chem.MolFromSmarts("[C@](O)(COP(=O)(O)OCCN)COC(=O)[C,H]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone with (R)-configuration and phosphate-ethanolamine linkage not found"

    # Ensure the 1-O-acyl group is linked in the correct position
    # This checks for an ester bond between acyl chain and glycerol
    acyl_pattern = Chem.MolFromSmarts("OC(=O)[C,H]")
    if not mol.HasSubstructMatch(acyl_pattern):
        return False, "No 1-O-acyl ester linkage found"
    
    # Exclude structural elements common to phospholipids like PC or PS
    # This should help differentiate PE from other phosphate-containing lipids
    wrong_substructs = [
        Chem.MolFromSmarts("O=C(OCC[N+](C)(C)C)P(=O)(O)O"), # Common pattern in PC
        Chem.MolFromSmarts("OC[C@H](NP(=O)(O)O)CO") # Pattern typical of PS
    ]
    for sub in wrong_substructs:
        if mol.HasSubstructMatch(sub):
            return False, "Found structural elements common to other phospholipids like PC or PS"
    
    return True, "Molecule fits all structural criteria for 1-acyl-sn-glycero-3-phosphoethanolamine"