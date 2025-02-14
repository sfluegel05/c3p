"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    Uses a SMARTS-based approach for ring detection and alpha configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for steroid core (generalized tetracyclic ring system)
    steroid_core_pattern = Chem.MolFromSmarts("[R1]1[R1][R1][R1]2[R1][R1][R1]3[R1]([R1][R1]4[R1][R1]1[R1]2)[R1][R1]43")
    if steroid_core_pattern is None:
        return None, "Invalid steroid core SMARTS pattern"
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Not a steroid core structure"

    # 2. Check for a 3-alpha-hydroxy group with chirality
    alpha_hydroxy_pattern = Chem.MolFromSmarts("[C@H](O)[C]")
    if alpha_hydroxy_pattern is None:
        return None, "Invalid alpha hydroxy SMARTS pattern"
    
    alpha_matches = mol.GetSubstructMatches(alpha_hydroxy_pattern)
    
    if len(alpha_matches) == 0:
        
        #Check if it is a beta hydroxy group for a more descriptive negative outcome.
        beta_hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)[C]")
        if beta_hydroxy_pattern is None:
            return None, "Invalid beta hydroxy SMARTS pattern"

        beta_matches = mol.GetSubstructMatches(beta_hydroxy_pattern)
        if len(beta_matches) > 0:
            return False, "Hydroxyl group at position 3 is in beta-configuration"

        return False, "No hydroxyl group found at position 3"

    # Check that the hydroxyl group is attached to the steroid core by matching the pattern with the molecule
    matched_atoms = mol.GetSubstructMatch(alpha_hydroxy_pattern)
    if not matched_atoms:
        return False, "No hydroxyl group attached to steroid core"

    core_match = mol.GetSubstructMatch(steroid_core_pattern)
    
    if not core_match:
         return False, "No steroid core match after hydroxy match. Internal Error"

    hydroxyl_carbon_index = matched_atoms[0]
    
    is_core_attached = False
    for core_atom_index in core_match:
        
        #Get all neighbors of the core atoms
        for neighbor_atom in mol.GetAtomWithIdx(core_atom_index).GetNeighbors():
            if neighbor_atom.GetIdx() == hydroxyl_carbon_index:
                is_core_attached = True
                break
        if is_core_attached:
            break


    if not is_core_attached:
        return False, "Hydroxyl group not attached to the steroid core"

    return True, "Molecule is a 3alpha-hydroxy steroid"