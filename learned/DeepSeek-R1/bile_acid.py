"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are hydroxy-5beta-cholanic acids with a steroid nucleus and specific stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal carboxylic acid (position 24)
    # Pattern: 4 carbons in chain ending with COOH (C-24)
    side_chain_acid = Chem.MolFromSmarts("[CX4]-[CX4]-[CX4]-[CX4]-C(=O)O")
    if not mol.HasSubstructMatch(side_chain_acid):
        return False, "Missing 24-oic acid side chain"

    # Verify steroid nucleus (gonane core: three 6-membered rings + one 5-membered)
    steroid_core = Chem.MolFromSmarts(
        "[C@H]1[C@@H]2CC[C@H]3[C@]4(CCC5CC(C(=O)O)CC[C@]5(C)[C@H]4CC[C@]3([C@H]2C[C@@H]1C)C)C"
    )
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid nucleus detected"

    # Check 5β configuration (trans A/B ring fusion)
    # Pattern for 5β: C5 has R configuration (axial H)
    ab_fusion = Chem.MolFromSmarts("[C@H]1[C@@H](C)[C@@H]2CC[C@H]3[C@@]")
    if not mol.HasSubstructMatch(ab_fusion):
        return False, "5beta configuration not confirmed"

    # Check for at least one hydroxyl on steroid nucleus (positions 3,7,12 etc.)
    hydroxyl_pattern = Chem.MolFromSmarts("[C][OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if not hydroxyl_matches:
        return False, "No hydroxyl groups on steroid core"

    # Verify side chain attachment at C17 (D-ring)
    d_ring_attachment = Chem.MolFromSmarts("[C@H]1CC[C@H]([C@H](CCCC(=O)O)C)[C@@]1")
    if not mol.HasSubstructMatch(d_ring_attachment):
        return False, "Side chain not attached at C17"

    return True, "5β-cholanic acid with steroid nucleus, hydroxyl groups, and 24-oic acid"