"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    Must have the myo-inositol core structure (specific stereochemistry) with at least one phosphate group.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Define myo-inositol core SMARTS with correct stereochemistry (axial OH at position 2)
    # SMARTS based on myo-inositol structure: O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O
    myo_core_smarts = Chem.MolFromSmarts("""
        [C@H]1([C@H]([C@@H]([C@H]([C@@H]([C@@H]1O)O)O)O)O)O
    """)
    if not mol.HasSubstructMatch(myo_core_smarts):
        return False, "No myo-inositol core found"

    # Find all oxygen atoms attached to the core (potential hydroxyl/phosphate sites)
    core_match = mol.GetSubstructMatch(myo_core_smarts)
    core_o = []
    for atom_idx in core_match:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetIdx() not in core_match:
                core_o.append(neighbor.GetIdx())

    # Check for at least one phosphate group attached to core oxygen
    phosphate_pattern = Chem.MolFromSmarts("[O][P](=O)([O])[O]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate groups detected"

    # Verify at least one phosphate is attached to the inositol core
    core_attached_phosphate = False
    for match in phosphate_matches:
        p_idx = match[0]
        # Check if this phosphate is connected to a core oxygen
        for o_idx in core_o:
            if mol.GetBondBetweenAtoms(o_idx, p_idx) is not None:
                core_attached_phosphate = True
                break
        if core_attached_phosphate:
            break

    if not core_attached_phosphate:
        return False, "Phosphate not attached to myo-inositol core"

    # Check all core substituents are either hydroxyl or phosphate groups
    valid_substituents = True
    for o_idx in core_o:
        atom = mol.GetAtomWithIdx(o_idx)
        # Check what's attached to the oxygen (should be H or P)
        if atom.GetDegree() == 1:  # Hydroxyl (-OH)
            continue
        # Check if connected to phosphorus (phosphate)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 15:
                break
        else:
            valid_substituents = False
            break

    if not valid_substituents:
        return False, "Non-phosphate/hydroxyl substituents present"

    return True, "myo-inositol core with phosphate group(s)"