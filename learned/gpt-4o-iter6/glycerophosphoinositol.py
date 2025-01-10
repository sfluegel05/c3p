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
    
    # Define patterns more inclusively and logically
    glycerol_pattern = Chem.MolFromSmarts("[O][C@H](CO)C([O])")  # Assume correct stereochemistry from input
    inositol_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")  # More explicit phosphate pattern

    # Check for the presence of glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Search for inositol ring and ensure connection to a phosphate group at the sn-3 position
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    phosphate_connected = False
    
    for match in inositol_matches:
        # For each match, check for nearby phosphate connection
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'P': # Look for a phosphate directly bonded to inositol
                phosphate_connected = True
                break
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'P':
                    phosphate_connected = True
                    break
            if phosphate_connected:
                break

    if not phosphate_connected:
        return False, "Phosphate group is not attached to inositol at the correct position"
    
    # Check for at least one ester linkage to fatty acids
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]CO")  # Assume correct stereochemistry from input
    if len(mol.GetSubstructMatches(ester_pattern)) < 1:
        return False, "No ester linkages found, expected fatty acids"

    return True, "Contains glycerophosphoinositol structure (glycerol backbone, attached phosphate group to inositol, inositol, and fatty acids)"