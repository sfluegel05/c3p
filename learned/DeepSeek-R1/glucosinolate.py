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
    1) Thioglucose moiety (sulfur attached to glucose)
    2) Central carbon connected to sulfur (from thioglucose) and sulfonated oxime group (C=N-OSO3^-)
    3) Anti configuration between side chain and sulfate group
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # 1. Check thioglucose pattern (sulfur connected to glucose backbone)
    thiogluc_pattern = Chem.MolFromSmarts("[C@H]1(O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)CO)S")
    if not mol.HasSubstructMatch(thiogluc_pattern):
        return False, "Missing thioglucose moiety"

    # 2. Find sulfonated oxime group (C=N-OSO3^-)
    # Pattern accounts for different protonation states and resonance
    oxime_sulf_pattern = Chem.MolFromSmarts("[CX3]=[NX2]-[OX2]-S(=O)(=O)[O-]")
    oxime_matches = mol.GetSubstructMatches(oxime_sulf_pattern)
    if not oxime_matches:
        return False, "Missing sulfonated oxime group"

    # 3. Verify sulfur from thioglucose is bonded to the central carbon (C in C=N-O-SO3)
    # Get indices: [C]=N-O-SO3^-
    central_c = [match[0] for match in oxime_matches]
    thiogluc_s = [a.GetBeginAtomIdx() for a in mol.GetSubstructMatches(thiogluc_pattern) if mol.GetAtomWithIdx(a[0]).GetAtomicNum() == 16]

    # Check connectivity between thioglucose S and central C
    connected = False
    for s_idx in thiogluc_s:
        s_atom = mol.GetAtomWithIdx(s_idx)
        for neighbor in s_atom.GetNeighbors():
            if neighbor.GetIdx() in central_c:
                connected = True
                break
        if connected:
            break
    if not connected:
        return False, "Thioglucose S not bonded to central C=N carbon"

    # 4. Check anti configuration using SMARTS pattern (approximate)
    # Looks for / and \ in the C=N-O-SO3 and adjacent chain
    anti_pattern = Chem.MolFromSmarts("[SX2]/C(=N/O/S(=O)(=O)[O-])/C")
    if not mol.HasSubstructMatch(anti_pattern):
        anti_pattern2 = Chem.MolFromSmarts("[SX2]\\C(=N\\O/S(=O)(=O)[O-])/C")
        if not mol.HasSubstructMatch(anti_pattern2):
            return False, "Anti configuration not detected"

    # 5. Check for side chain (R group) on central carbon
    # At least one non-oxygen/nitrogen/sulfur atom attached to central C
    r_group = False
    for c_idx in central_c:
        c_atom = mol.GetAtomWithIdx(c_idx)
        for neighbor in c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() not in {7,8,16}:
                r_group = True
                break
        if r_group:
            break
    if not r_group:
        return False, "Missing side chain on central carbon"

    return True, "Contains thioglucose, sulfonated oxime, correct connectivity and stereochemistry"