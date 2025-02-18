"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: CHEBI:16204 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    Must have: myo-inositol core with specific stereochemistry, phosphate at position 1,
    and two fatty acid chains attached via sn-glycero-3-phospho linkage.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for myo-inositol core with correct stereochemistry
    # Pattern matches 1D-myo-inositol with phosphate at position 1
    inositol_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([C@H]([C@@H]([C@H]([C@H]1O)O)O)O)O)OP(=O)(O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "Missing or incorrect myo-inositol core"

    # Verify phosphate is connected to glycerol backbone
    phosphate_glycerol_pattern = Chem.MolFromSmarts("[C@H](COP(=O)(O)O[C@])(COC(=O))COC(=O)")
    if not mol.HasSubstructMatch(phosphate_glycerol_pattern):
        return False, "Incorrect phosphate-glycerol linkage"

    # Check for two esterified fatty acid chains (at least 4 carbons each)
    ester_pattern = Chem.MolFromSmarts("[CX4]-[OX2]-C(=O)-[CX4]")
    esters = mol.GetSubstructMatches(ester_pattern)
    if len(esters) < 2:
        return False, f"Found {len(esters)} ester groups, need at least 2"

    # Verify chain lengths (crude check for at least 4 carbons in each chain)
    for ester in esters:
        chain_start = ester[2]  # Carbon in carbonyl group
        chain_length = 0
        atom = mol.GetAtomWithIdx(chain_start)
        while True:
            neighbors = [n for n in atom.GetNeighbors() 
                        if n.GetIdx() != ester[1] and n.GetSymbol() == "C"]
            if not neighbors:
                break
            atom = neighbors[0]
            chain_length += 1
        if chain_length < 3:  # 3 carbons + carbonyl = 4 total
            return False, "Fatty acid chain too short"

    # Check molecular weight (typical >700 Da)
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for phosphatidylinositol"

    return True, "Contains 1D-myo-inositol with phosphatidyl group at position 1"