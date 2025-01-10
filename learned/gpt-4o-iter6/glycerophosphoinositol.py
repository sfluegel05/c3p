"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem

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
    
    # Enhanced glycerol backbone pattern (include stereochemistry possibilities)
    glycerol_pattern = Chem.MolFromSmarts("O[C@@H](CO)C(O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Define inositol and phosphate patterns
    inositol_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)")
    
    # Search for inositol ring connected to a phosphate group
    inositol_matches = mol.HasSubstructMatch(inositol_pattern)
    if not inositol_matches:
        return False, "Inositol ring not found"
    
    phosphate_matches = mol.HasSubstructMatch(phosphate_pattern)
    if not phosphate_matches:
        return False, "Phosphate group not found"

    # Ensure phosphate is connected to inositol at sn-3 position
    phosphate_connected = False
    phosphate_atoms = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'P'}
    for match in mol.GetSubstructMatches(inositol_pattern):
        if phosphate_atoms.intersection(match):
            phosphate_connected = True
            break
    
    if not phosphate_connected:
        return False, "Phosphate group is not attached to inositol"
    
    # Check for at least one ester linkage to indicate fatty acids
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]CO")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester linkages found, expected fatty acids"

    return True, "Contains glycerophosphoinositol structure (glycerol backbone, attached phosphate group at sn-3, inositol, and fatty acids)"