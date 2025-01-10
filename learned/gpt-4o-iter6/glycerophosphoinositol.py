"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol is a glycerophospholipid having the polar alcohol inositol esterified
    to the phosphate group at the sn-3 position of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Glycerol backbone substructure (simple form)
    glycerol_pattern = Chem.MolFromSmarts("[OH][CH2][CH]([OH])[CH2][OH]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Phosphate group esterified with inositol
    phosphate_inositol_pattern = Chem.MolFromSmarts("P(=O)(O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)(OC)=O") # Simplified inositol representation
    if not mol.HasSubstructMatch(phosphate_inositol_pattern):
        return False, "Phosphate inositol unit in the required format not found"
    
    # Check for fatty acid ester linkages 
    ester_pattern = Chem.MolFromSmarts("O=C[O][CH2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "Less than 1 ester linkage found"
    
    # Check that phosphate group with inositol is at sn-3 position
    inositol_sn3_pattern = Chem.MolFromSmarts("[CH2][O][P](=O)(O)[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(inositol_sn3_pattern):
        return False, "Phosphate inositol not attached at sn-3 position"

    return True, "Contains a glycerol backbone with a phosphate group and inositol esterified at sn-3 position"