"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component
    has myo-configuration and at least one of the hydroxyl groups is replaced with a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        (bool, str): True if the molecule is a myo-inositol phosphate with the reason, else False with reason.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for myo-inositol core
    myo_inositol_pattern = Chem.MolFromSmarts("C1(O)C(O)C(O)C(O)C(O)C1")
    
    # Check for presence of myo-inositol core
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol core structure found"
    
    # Define SMARTS pattern for phosphate group
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    
    # Get matches for phosphate groups and count them
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate groups found"
    
    # Ensure direct attachment of phosphate to the inositol core
    phosphate_attached = False
    for match in phosphate_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'O':  # We assume P=O connection.
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:  # Check if P-O connects to inositol's carbon.
                        cycle_info = neighbor.IsInRingSize(6)
                        if cycle_info:  # Check for membership in a 6-membered ring, characteristic of inositol
                            phosphate_attached = True
                            break
                if phosphate_attached:
                    break
        if phosphate_attached:
            break
            
    if not phosphate_attached:
        return False, "Phosphate not directly attached to inositol core"

    return True, "Valid myo-inositol phosphate with correctly placed phosphate groups"