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

    # Define myo-inositol core SMARTS with correct stereochemistry
    # SMARTS matches six-membered ring with specific chiral centers and oxygen substituents
    myo_core_smarts = Chem.MolFromSmarts("""
        [C@H]1([OX2])
        [C@H]([OX2])
        [C@@H]([OX2])
        [C@H]([OX2])
        [C@@H]([OX2])
        [C@@H]1[OX2]
    """)
    if myo_core_smarts is None:
        return False, "Failed to parse myo-inositol SMARTS"
    
    if not mol.HasSubstructMatch(myo_core_smarts):
        return False, "No myo-inositol core found"

    # Get core oxygen atoms (attached to the 6 carbons)
    core_match = mol.GetSubstructMatch(myo_core_smarts)
    core_carbons = set(core_match)
    core_oxygens = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetDegree() > 0:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetIdx() in core_carbons:
                core_oxygens.append(atom.GetIdx())

    # Check for at least one phosphate group attached to core oxygen
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)[OX2]")
    if phosphate_pattern is None:
        return False, "Failed to parse phosphate SMARTS"
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    
    has_phosphate = False
    for match in phosphate_matches:
        # Check if any oxygen in the phosphate is a core oxygen
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8 and atom.GetIdx() in core_oxygens:
                has_phosphate = True
                break
        if has_phosphate:
            break
    
    if not has_phosphate:
        return False, "No phosphate attached to myo-inositol core"

    # Verify all core substituents are hydroxyl or phosphate
    for o_idx in core_oxygens:
        atom = mol.GetAtomWithIdx(o_idx)
        if atom.GetDegree() == 1:  # Hydroxyl group
            continue
        # Check if connected to phosphorus in a phosphate group
        has_phosphorus = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 15:
                has_phosphorus = True
                break
        if not has_phosphorus:
            return False, "Non-phosphate substituent on core"

    return True, "myo-inositol core with phosphate group(s)"