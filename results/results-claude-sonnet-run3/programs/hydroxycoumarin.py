from rdkit import Chem
from rdkit.Chem import AllChem

def is_hydroxycoumarin(smiles: str):
    """
    Determines if a molecule is a hydroxycoumarin (coumarin with at least one hydroxy substituent).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a hydroxycoumarin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First find the coumarin core
    # Using a more general SMARTS pattern for 2H-chromen-2-one (coumarin) core
    coumarin_pattern = Chem.MolFromSmarts('c12ccccc1OC(=O)C=C2')
    
    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin core structure found"

    # Find all hydroxyl groups
    hydroxy_pattern = Chem.MolFromSmarts('[OH]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    if not hydroxy_matches:
        return False, "No hydroxyl groups found"

    # Get the atoms that are part of the coumarin core
    coumarin_atoms = set(sum(mol.GetSubstructMatches(coumarin_pattern), ()))

    # Check if any hydroxyl group is attached to the coumarin system
    hydroxy_on_coumarin = False
    for match in hydroxy_matches:
        oh_atom = match[0]
        # Get the atom the OH is connected to
        for neighbor in mol.GetAtomWithIdx(oh_atom).GetNeighbors():
            if neighbor.GetIdx() in coumarin_atoms:
                hydroxy_on_coumarin = True
                break
        if hydroxy_on_coumarin:
            break

    if not hydroxy_on_coumarin:
        return False, "Hydroxyl groups present but not on coumarin core"

    # Count total number of hydroxyl groups
    total_hydroxy = len(hydroxy_matches)
    
    return True, f"Found coumarin core with {total_hydroxy} hydroxyl group(s)"
# Pr=None
# Recall=0.0