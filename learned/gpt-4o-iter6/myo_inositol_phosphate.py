"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component
    has myo-configuration and at least one hydroxyl group is replaced with a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        (bool, str): True if the molecule is a myo-inositol phosphate with the reason, else False with reason.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for myo-inositol core with correct stereochemistry
    myo_inositol_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@@H](O)[C@H]1O")
    
    # Check for presence of an inositol core
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol core structure with correct stereochemistry found"
    
    # Define SMARTS pattern for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    
    # Get matches for phosphate groups
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate groups found"

    # Ensure that at least one phosphate is directly connected to the inositol core
    for phosphate_match in phosphate_matches:
        phosphate_atom = mol.GetAtomWithIdx(phosphate_match[0])
        # Ensure phosphate group is directly attached to a carbon in the inositol ring
        attached_to_inositol = any(
            neighbor.IsInRingSize(6) and neighbor.GetAtomicNum() == 6 
            for neighbor in phosphate_atom.GetNeighbors()
        )
        if attached_to_inositol:
            return True, "Valid myo-inositol phosphate with correctly attached phosphate groups"

    return False, "Phosphate not directly attached to inositol core"