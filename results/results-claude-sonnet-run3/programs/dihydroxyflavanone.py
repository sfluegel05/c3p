from rdkit import Chem
from rdkit.Chem import AllChem

def is_dihydroxyflavanone(smiles: str):
    """
    Determines if a molecule is a dihydroxyflavanone (flavanone with exactly 2 hydroxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroxyflavanone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for flavanone core structure
    # SMARTS pattern for flavanone core - more flexible pattern
    flavanone_pattern = Chem.MolFromSmarts('O1CC(=O)c2c(cc[cR1]c2)C1')
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Not a flavanone - missing core structure"

    # Count number of hydroxy groups
    oh_pattern = Chem.MolFromSmarts('[OH]')
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    # Count only hydroxy groups connected to aromatic carbons
    true_oh_count = 0
    for match in oh_matches:
        oh_atom = mol.GetAtomWithIdx(match[0])
        neighbors = oh_atom.GetNeighbors()
        if len(neighbors) == 1:  # Make sure it's terminal OH
            neighbor = neighbors[0]
            # Check if the neighbor is aromatic carbon
            if neighbor.GetIsAromatic() and neighbor.GetSymbol() == 'C':
                true_oh_count += 1

    if true_oh_count != 2:
        return False, f"Has {true_oh_count} hydroxy groups instead of required 2"

    # If we get here, it's a flavanone with exactly 2 hydroxy groups
    return True, "Contains flavanone core with exactly 2 hydroxy groups"
# Pr=None
# Recall=0.0