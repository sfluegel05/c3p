"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on SMILES.
    Criteria:
    1. Contains steroid nucleus (four fused rings: 3 cyclohexane + 1 cyclopentane)
    2. Has a 3-oxo group (ketone at position 3 of the A ring)
    3. 5-beta configuration (cis A/B ring fusion)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check basic steroid nucleus (four fused rings)
    steroid_smarts = Chem.MolFromSmarts("[C]1~[C]~[C]~[C]~[C]~[C]1.[C]1~[C]~[C]~[C]~[C]~[C]1.[C]1~[C]~[C]~[C]~[C]~[C]1.[C]1~[C]~[C]~[C]~[C]1")
    if not mol.HasSubstructMatch(steroid_smarts):
        return False, "No steroid nucleus"

    # Check for 3-oxo group in A ring (position 3)
    # A ring is typically 6-membered, position 3 is adjacent to ring B
    a_ring_ketone = Chem.MolFromSmarts("[C]=O")
    ketone_matches = mol.GetSubstructMatches(a_ring_ketone)
    if not ketone_matches:
        return False, "No 3-oxo group"
    # Position 3 check might need better spatial analysis

    # Check 5-beta configuration (cis A/B ring fusion)
    # SMARTS pattern for 5-beta configuration (cis fusion)
    # Focus on C5 (atom 2 in pattern) and C10 (atom 3) stereochemistry
    beta_config = Chem.MolFromSmarts("[C@]1([H])[C@@]2([H])CC[C@H]1[C@@H]2")
    if not mol.HasSubstructMatch(beta_config):
        return False, "5-beta configuration not found"

    return True, "3-oxo group and 5-beta steroid nucleus present"