"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is typically acknowledged as a sugar alcohol with multiple hydroxyl groups and no reactive carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecular structure
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the presence of carbonyl groups
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")  # C=O in aldehydes and ketones
    if mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Contains carbonyl-like group(s), disqualifying it as a pure alditol"

    # Count hydroxyl groups (attached as -OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, f"Found {len(hydroxyl_matches)} hydroxyl groups, need at least 3 for a typical alditol"

    # Check for conformity of alditol-specific traits (simple hydrocarbon chains with OHs only)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8, 1):  # Carbon, Oxygen, Hydrogen
            return False, f"Contains atoms other than C, O, and H: {atom.GetSymbol()}"
    
    # Validation allowing for some ring structures (as long as the core is alditol-like)
    alditol_like = True
    for ring in mol.GetRingInfo().AtomRings(): 
        # Criteria could be the presence of glycosidic-like linkages or beyond simple cyclic alditol
        if len(ring) > 6 or not any(mol.GetAtomWithIdx(at).GetSymbol() == 'O' for at in ring):
            alditol_like = False
            break

    if not alditol_like:
        return False, "Structure contains complex rings that abstract from a simple alditol classification"

    return True, "Contains multiple hydroxyl groups and appropriate structure, consistent with an alditol-like compound"