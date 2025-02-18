"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: CHEBI:77764 11beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    Requirements: Must have steroid nucleus and beta-oriented hydroxyl at position 11.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid check - look for characteristic fused ring system
    # This SMARTS pattern matches the core tetracyclic structure of steroids
    steroid_pattern = Chem.MolFromSmarts("[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Missing steroid nucleus"

    # 11-beta-OH check - look for hydroxyl in specific stereochemical environment
    # Pattern matches hydroxyl at position 11 with beta orientation (same face as C18 methyl)
    beta_oh_pattern = Chem.MolFromSmarts("[C@@H]([OH])[C@H]1CC[C@H]2[C@@H]3CC(=O)CC[C@H]3C[C@H]21")
    if not mol.HasSubstructMatch(beta_oh_pattern):
        return False, "No 11-beta-hydroxyl group detected"

    return True, "Contains steroid nucleus with 11-beta-hydroxyl group"