"""
Classifies: CHEBI:24279 glucosinolate
"""
"""
Classifies: CHEBI:24163 glucosinolate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    Must have:
    1) Thioglucose moiety (sulfur attached to glucose backbone)
    2) Central carbon connected to sulfur (from thioglucose) and sulfonated oxime group (C=N-O-SO3^-)
    3) Anti configuration between side chain and sulfonate group
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # 1. Check for thioglucose moiety (S connected to glucose-like structure)
    # Generalized pattern for sulfur attached to hexose (C1O... with S connected)
    thiogluc_pattern = Chem.MolFromSmarts("[C][C@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)CO")
    thiogluc_matches = mol.GetSubstructMatches(thiogluc_pattern)
    if not thiogluc_matches:
        return False, "Missing thioglucose moiety"

    # Get sulfur atoms attached to the glucose moiety
    thiogluc_s = []
    for match in thiogluc_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 16:  # Sulfur
                thiogluc_s.append(atom_idx)
    if not thiogluc_s:
        return False, "No sulfur in thioglucose moiety"

    # 2. Find sulfonated oxime group (C=N-O-SO3^-)
    oxime_sulf_pattern = Chem.MolFromSmarts("[CX3]=[NX2]-[OX2]-S(=O)(=O)[O-]")
    oxime_matches = mol.GetSubstructMatches(oxime_sulf_pattern)
    if not oxime_matches:
        # Try alternative pattern with possible different charge states
        oxime_sulf_pattern = Chem.MolFromSmarts("[CX3]=[NX2]-[OX2]S(=O)(=O)[O-]")
        oxime_matches = mol.GetSubstructMatches(oxime_sulf_pattern)
        if not oxime_matches:
            return False, "Missing sulfonated oxime group"

    # 3. Verify sulfur from thioglucose is bonded to the central C of the oxime group
    connected = False
    for (c_idx, n_idx, o_idx, s_idx) in oxime_matches:
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Check if this C is bonded to the thioglucose S
        for neighbor in c_atom.GetNeighbors():
            if neighbor.GetIdx() in thiogluc_s:
                connected = True
                break
        if connected:
            break
    if not connected:
        return False, "Thioglucose S not bonded to central C=N carbon"

    # 4. Check stereochemistry of C=N bond (must have double bond stereochemistry for anti)
    # Look for either / or \ in the C=N-O-SO3 group
    stereo_bond_present = False
    for (c_idx, n_idx, o_idx, s_idx) in oxime_matches:
        bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
        if bond.GetStereo() != Chem.BondStereo.STEREONONE:
            stereo_bond_present = True
            break
    if not stereo_bond_present:
        return False, "Missing stereochemistry on C=N bond (anti config required)"

    # 5. Check for side chain (R group) attached to central C
    has_side_chain = False
    for (c_idx, n_idx, o_idx, s_idx) in oxime_matches:
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Should have at least two bonds: one to S, one to R group
        if len(c_atom.GetBonds()) >= 2:
            for bond in c_atom.GetBonds():
                neighbor = bond.GetOtherAtom(c_atom)
                if neighbor.GetAtomicNum() not in {7, 8, 16}:  # Not N, O, S
                    has_side_chain = True
                    break
            if has_side_chain:
                break
    if not has_side_chain:
        return False, "Missing side chain on central carbon"

    return True, "Contains thioglucose, sulfonated oxime, correct connectivity and stereochemistry"