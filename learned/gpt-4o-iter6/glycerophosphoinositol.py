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
    
    # Define patterns for key structural components
    glycerol_pattern = Chem.MolFromSmarts("O[C@H](COP)C(O)")  # Glycerol backbone with phosphate
    inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1O")
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")

    # Check for the presence of glycerol backbone with phosphate
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with phosphate group found"

    # Check for inositol attached to phosphate
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No inositol group found"

    # Confirm phosphate linkage to inositol
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    
    for p_match in phosphate_matches:
        for i_match in inositol_matches:
            phosphate_atom = mol.GetAtomWithIdx(p_match[0])  # Assuming first atom of match is P
            if phosphate_atom.GetSymbol() == 'P' and \
               any(neighbor.GetIdx() in i_match for neighbor in phosphate_atom.GetNeighbors()):
                return True, "Contains glycerophosphoinositol structure (glycerol backbone, attached phosphate group to inositol, inositol, and fatty acids)"

    return False, "Phosphate group is not correctly connected to inositol"