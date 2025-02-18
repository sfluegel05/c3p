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

    # Check for 1D-myo-inositol core with correct stereochemistry and phosphate at position 1
    # Inositol pattern: six-membered ring with specific hydroxyl configurations and phosphate at O1
    inositol_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]([C@H]([C@@H]([C@H]([C@H]1O)O)O)O)OP(=O)(O)O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "Missing or incorrect 1D-myo-inositol core with phosphate at position 1"

    # Verify sn-glycero-3-phospho linkage: glycerol with two esters (sn-1 and sn-2) and phosphate at sn-3
    # Glycerol pattern: central carbon (sn-2) has R configuration ([C@@H]) with two esters and phosphate
    glycerol_pattern = Chem.MolFromSmarts("[C@@H](-O-C(=O))(-[CH2]-O-C(=O))-[CH2]-O-P(=O)(O)-O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Incorrect sn-glycero-3-phospho linkage to inositol"

    # Check for exactly two fatty acid ester groups (O-C=O connected to carbon chains)
    ester_matches = mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2]C(=O)"))
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Remove chain length check to accommodate short chains (e.g., dibutyrl example)
    # Check molecular weight (adjust threshold as needed)
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 400:  # Lowered threshold to include smaller examples
        return False, f"Molecular weight too low ({mol_wt:.1f} Da)"

    return True, "Contains 1D-myo-inositol with sn-glycero-3-phospho and two fatty acid chains"