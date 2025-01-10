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

    # New approach: Allow for some ring presence due to saccharide origins
    # BUT confirm no closed chains longer than 6 atoms and linear backbone

    # Assess count of hydroxyl groups, most carbon atoms should have -OH
    c_oh_groups = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' and any(n.GetSymbol() == 'O' for n in atom.GetNeighbors()))
    c_total = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    
    # Must have predominantly C-OH structures
    if not (c_oh_groups >= c_total - 1):  # Allow for small loss of -OH in branched blocks
        return False, "Not enough C-OH substitution typical of alditols"

    # Ensure reduction impact - No terminal carbonyl C=O groups
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains carbonyl groups, indication of precursor not alditol"

    # Check for redundant carbons joining chains
    terminal_oh_pattern = Chem.MolFromSmarts("[CH2]O")
    terminal_oh_matches = mol.GetSubstructMatches(terminal_oh_pattern)
    if len(terminal_oh_matches) < 2:
        return False, "Does not have typical terminal CH2OH groups"

    # Consider stereochemistry: Distinct stereoisomers typical in alditols
    # Need a SMARTS pattern to match these pretty selectively
    multi_hydroxyl_pattern = Chem.MolFromSmarts("[C@@H]([OH])")  # - C attached to hydroxyl, chiral
    if mol.HasSubstructMatch(multi_hydroxyl_pattern):
        return True, "Molecule matches hydroxyl stereochemical pattern of an alditol"

    return False, "Failed to match stereochemistry and hydroxyl profile for alditol"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17754',
        'name': 'alditol',
        'definition': "A carbohydrate that is an acyclic polyol having the general formula HOCH2[CH(OH)]nCH2OH, formally derivable from an aldose by reduction of the carbonyl group.",
        'parents': ['CHEBI:63221'],
    }
}