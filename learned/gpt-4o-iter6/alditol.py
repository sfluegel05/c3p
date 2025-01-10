"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula 
    HOCH2[CH(OH)]nCH2OH (formally derived from an aldose by reduction).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for acyclic structure: no rings present
    if mol.GetRingInfo().IsInitialized() and mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains cyclic structure"
    
    # Check for polyol structure: multiple C-OH groups
    hydroxyl_pattern = Chem.MolFromSmarts("[C;H2,H3][OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Not enough C-OH groups to be a polyol"
    
    # Check for terminal OH groups (CH2OH ends)
    terminal_oh_pattern = Chem.MolFromSmarts("[CX3](O)(H)H")
    terminal_oh_matches = mol.GetSubstructMatches(terminal_oh_pattern)
    if len(terminal_oh_matches) < 2:
        return False, "Does not have terminal CH2OH groups"
    
    # Verify absence of carbonyl groups (C=O) indicating reduction
    carbonyl_pattern = Chem.MolFromSmarts("[C]=[O]")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains carbonyl (C=O) group, not an alditol"
    
    return True, "Molecule matches structure of an alditol"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17754',
        'name': 'alditol',
        'definition': "A carbohydrate that is an acyclic polyol having the general formula HOCH2[CH(OH)]nCH2OH, formally derivable from an aldose by reduction of the carbonyl group.",
        'parents': ['CHEBI:63221'],
    }
}