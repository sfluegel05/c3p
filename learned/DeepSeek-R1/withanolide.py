"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: CHEBI:??? withanolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is a steroid lactone with a modified side chain forming a lactone ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for lactone ring (cyclic ester)
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[O]C(=O)'))
    lactone_found = False
    for match in ester_matches:
        o_idx, c_idx = match[0], match[1]
        # Check if O and C=O are in the same ring
        ri = mol.GetRingInfo()
        for ring in ri.AtomRings():
            if o_idx in ring and c_idx in ring:
                lactone_found = True
                break
        if lactone_found:
            break
    if not lactone_found:
        return False, "No lactone ring found"

    # Check for steroid nucleus (tetracyclic system)
    # SMARTS pattern for steroid nucleus (simplified)
    steroid_pattern = Chem.MolFromSmarts(
        "[C@H]1CC[C@H]2[C@@H]3CC[C@H]4[C@H](CC[C@@H]4)[C@@H]3CC[C@@H]12"
    )
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"

    # Check if lactone is part of the side chain (approximate check)
    # Get atoms in steroid nucleus
    steroid_atoms = set()
    for match in mol.GetSubstructMatches(steroid_pattern):
        steroid_atoms.update(match)
    # Check if any atom in the lactone ring is not part of the steroid nucleus
    # (assuming lactone is in the side chain)
    lactone_in_sidechain = False
    for ring in ri.AtomRings():
        if lactone_found and any(atom in ring for atom in steroid_atoms):
            # Check if any atom in the lactone ring is not in the steroid nucleus
            if any(atom not in steroid_atoms for atom in ring):
                lactone_in_sidechain = True
                break
    if not lactone_in_sidechain:
        return False, "Lactone not in side chain"

    return True, "Steroid nucleus with lactone ring in side chain"