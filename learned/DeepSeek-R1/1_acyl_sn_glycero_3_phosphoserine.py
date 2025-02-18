"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
"""
Classifies: CHEBI:17544 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    Must have: - sn-glycerol backbone with phosphoserine at position 3
               - One acyl ester at position 1
               - Free hydroxyl at position 2

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches structure, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Core structure pattern with stereochemistry requirements
    # [1] Acyl ester at position 1 (O-C(=O)-R)
    # [4] sn-2 carbon with hydroxyl (specific stereochemistry)
    # [5] Phosphate at position 3 connected to serine
    core_pattern = Chem.MolFromSmarts(
        "[O][C](=O)-[O]-[CH2][C@@H](O)[CH2]-O-P(=O)(-O-C[C@H](N)C(=O)O)-O"
    )
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core sn-glycerophosphoserine structure not found"

    # Verify only one acyl group (exactly one ester linkage to glycerol)
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[CH2]-O-C(=O)"))
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} acyl groups, need exactly 1 at position 1"

    # Check sn-2 hydroxyl exists and no other esters
    sn2_carbon = core_pattern.GetAtomWithIdx(4)  # Index from SMARTS pattern
    for match in mol.GetSubstructMatches(core_pattern):
        sn2_idx = match[4]
        sn2_atom = mol.GetAtomWithIdx(sn2_idx)
        oxygen_count = 0
        for neighbor in sn2_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(sn2_idx, neighbor.GetIdx()).GetBondType() == Chem.BondType.SINGLE:
                oxygen_count += 1
        if oxygen_count != 1:
            return False, "sn-2 position must have exactly one hydroxyl group"

    return True, "1-acyl group on sn-glycerol-3-phosphoserine with correct stereochemistry"