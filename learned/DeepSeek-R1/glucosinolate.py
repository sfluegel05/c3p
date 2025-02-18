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
    A glucosinolate has a thioglucose moiety linked via a sulfur atom to a sulfonated oxime group
    with an anti configuration between the side chain and sulfate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glucosinolate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for thioglucose moiety (sulfur attached to a glucose-like structure)
    # Generalized pattern for sulfur connected to a hexose ring with multiple hydroxyls and hydroxymethyl
    thioglucose_pattern = Chem.MolFromSmarts(
        "[C@H]1(O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)CO)S"
    )
    if not mol.HasSubstructMatch(thioglucose_pattern):
        return False, "No thioglucose moiety detected"

    # Check for sulfonated oxime group (N-O-SO3^-)
    sulf_oxime_pattern = Chem.MolFromSmarts("[CX3]=[NX2][OX2][S](=O)(=O)[O-]")
    sulf_oxime_matches = mol.GetSubstructMatches(sulf_oxime_pattern)
    if not sulf_oxime_matches:
        return False, "No sulfonated oxime group found"

    # Verify connectivity between thioglucose sulfur and the central carbon in C=N-O-SO3
    thioglucose_s = [match[0] for match in mol.GetSubstructMatches(thioglucose_pattern)]
    sulf_oxime_c = [match[0] for match in sulf_oxime_matches]  # CX3 in the pattern

    connected = False
    for s_idx in thioglucose_s:
        s_atom = mol.GetAtomWithIdx(s_idx)
        for neighbor in s_atom.GetNeighbors():
            if neighbor.GetIdx() in sulf_oxime_c:
                connected = True
                break
        if connected:
            break
    if not connected:
        return False, "Thioglucose sulfur not connected to C=N-O-SO3 group"

    # Check anti configuration (approximated via presence of directional bonds in SMILES)
    # This is a heuristic and may not cover all cases
    anti_pattern = Chem.MolFromSmarts("[S]/C(=N\\O/S(=O)(=O)[O-])/C")
    if not mol.HasSubstructMatch(anti_pattern):
        return False, "Anti configuration not detected between side chain and sulfate"

    # Check for a side chain (R group) attached to the central carbon
    # The central carbon is part of the sulf_oxime_matches' CX3
    # Ensure at least one non-S, non-O, non-N atom attached (crude check for R group)
    r_group_present = False
    for c_idx in sulf_oxime_c:
        c_atom = mol.GetAtomWithIdx(c_idx)
        for neighbor in c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() not in {8, 7, 16}:  # Exclude O, N, S
                r_group_present = True
                break
        if r_group_present:
            break
    if not r_group_present:
        return False, "No side chain (R group) attached to central carbon"

    return True, "Contains thioglucose, sulfonated oxime, correct connectivity, and anti configuration"