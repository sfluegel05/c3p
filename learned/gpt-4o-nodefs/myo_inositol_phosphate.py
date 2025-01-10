"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    The molecule must contain a myo-inositol core with phosphate groups directly attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Flexible SMARTS pattern for the myo-inositol core (cyclohexane ring with hydroxyl groups)
    inositol_pattern = Chem.MolFromSmarts("C1([C@@H](O)C(O)[C@@H](O)[C@H](O)C1O)O")
    
    # Check if the molecule contains the myo-inositol core
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol core structure found"

    # SMARTS pattern to represent mono-, di-, and tri-phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)O")
    
    # Check for presence of phosphate groups directly attached to the core
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups directly attached to inositol core"

    # Ensure phosphate groups are directly attached to each hydroxyl group
    for match in phosphate_matches:
        atom_indexes = set(match)
        hydroxyls = [atom.GetIdx() for atom in mol.GetAtoms() 
                     if atom.GetSymbol() == 'O' and 
                     len([nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() not in atom_indexes]) == 1]
        
        if not any(mol.GetBondBetweenAtoms(hydroxyl, hydroxyl + 1) for hydroxyl in hydroxyls):
            return False, "Phosphate groups not directly attached to inositol hydroxyl groups"
    
    return True, "Contains a myo-inositol core with phosphate groups directly attached"