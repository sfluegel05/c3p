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
    Must have: 1D-myo-inositol core with phosphate at position 1 connected to sn-glycero-3-phospho group,
    two fatty acid chains attached via ester bonds to the glycerol.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 1D-myo-inositol core with correct stereochemistry and phosphate at position 1
    # SMARTS pattern for 1D-myo-inositol-1-phosphate with specific hydroxyl configuration
    inositol_pattern = Chem.MolFromSmarts("[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1OP(=O)([OX2])[OX2]")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "Missing or incorrect 1D-myo-inositol core with phosphate at position 1"

    # Verify sn-glycero-3-phospho linkage with correct stereochemistry and two ester groups
    # SMARTS for glycerol part: [C@@H] (R configuration at sn-2) with two ester groups and phosphate at sn-3
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](-[OX2]-C(=O)-[!O])-[OX2]-C(=O)-[!O]-[CH2]-[OX2]-P(=O)([OX2])-[OX2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Incorrect sn-glycero-3-phospho linkage or ester groups"

    # Check that exactly two ester groups are attached to the glycerol's oxygens
    # Find all ester groups (O-C=O) connected to the glycerol's sn-1 and sn-2 positions
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2]-C(=O)"))
    glycerol_esters = 0
    for match in ester_matches:
        o_idx, c_idx = match
        # Check if the oxygen is part of the glycerol's ester groups
        for bond in mol.GetAtomWithIdx(o_idx).GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                neighbor = bond.GetOtherAtom(mol.GetAtomWithIdx(o_idx))
                if neighbor.GetAtomicNum() == 6 and neighbor.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    glycerol_esters += 1
                    break
    if glycerol_esters != 2:
        return False, f"Found {glycerol_esters} ester groups on glycerol, need exactly 2"

    return True, "Contains 1D-myo-inositol with sn-glycero-3-phospho and two fatty acid chains"