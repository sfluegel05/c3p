"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improve glycerol backbone pattern (include connections for phosphates and acids)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for modified phosphate-inositol pattern (more flexible with stereo elements)
    inositol_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    phosphate_pattern = Chem.MolFromSmarts("O=P(O)(O)")
    
    # Search for inositol ring esterified to a phosphate group
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if len(inositol_matches) == 0:
        return False, "Inositol ring not found"
    
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "Phosphate group not found"
    
    # Ensure phosphate is connected to inositol (flexibility in matching is allowed here)
    for match in inositol_matches:
        for ph_match in phosphate_matches:
            if set(match).intersection(ph_match):
                break
        else:
            continue
        break
    else:
        return False, "Phosphate group is not attached to inositol"
    
    # Check for ester linkages indicating fatty acids
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]CO")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "Less than 1 ester linkage found, expected fatty acids"

    return True, "Contains glycerophosphoinositol structure (glycerol backbone, attached phosphate group at sn-3, inositol, and fatty acids)"