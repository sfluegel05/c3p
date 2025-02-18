"""
Classifies: CHEBI:16749 1-phosphatidyl-1D-myo-inositol
"""
"""
Classifies: CHEBI:16204 1-phosphatidyl-1D-myo-inositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_1_phosphatidyl_1D_myo_inositol(smiles: str):
    """
    Determines if a molecule is a 1-phosphatidyl-1D-myo-inositol based on its SMILES string.
    Must have: 1D-myo-inositol core with phosphate at position 1 connected to sn-glycero-3-phospho group,
    two fatty acid chains attached via ester bonds to the glycerol.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 1D-myo-inositol core with correct stereochemistry
    # Pattern matches inositol with phosphate at position 1 (O1)
    inositol_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]([C@H]([C@@H]([C@H]([C@H]1O)O)O)O)OP(=O)(O)")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "Missing or incorrect 1D-myo-inositol core with phosphate at position 1"

    # Verify glycerol-phosphate linkage (sn-glycero-3-phospho)
    # Pattern matches glycerol with two ester groups and phosphate connected to inositol
    # Central carbon (sn-2) has @@ stereochemistry (R configuration)
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](COC(=O))(COP(=O)(O)O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)COC(=O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Incorrect sn-glycero-3-phospho linkage to inositol"

    # Check for two fatty acid ester groups (minimum 4 carbons each)
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4][OX2]C(=O)"))
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Basic chain length check (each chain >= 4 carbons including carbonyl)
    for ester in ester_matches:
        chain_start = ester[1]  # Oxygen atom in ester
        atom = mol.GetAtomWithIdx(chain_start)
        neighbor = [n for n in atom.GetNeighbors() if n.GetIdx() != ester[0]][0]  # Carbonyl carbon
        chain_length = 0
        current_atom = neighbor
        while True:
            next_carbons = [n for n in current_atom.GetNeighbors() 
                           if n.GetAtomicNum() == 6 and n.GetIdx() != current_atom.GetIdx()]
            if not next_carbons:
                break
            current_atom = next_carbons[0]
            chain_length += 1
        if chain_length < 3:  # 3 carbons + carbonyl = 4 total
            return False, "Fatty acid chain too short (minimum 4 carbons)"

    # Molecular weight check (typical >700 Da)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 600:  # Allow some flexibility for short chains in test cases
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"

    return True, "Contains 1D-myo-inositol with sn-glycero-3-phospho and two fatty acid chains"