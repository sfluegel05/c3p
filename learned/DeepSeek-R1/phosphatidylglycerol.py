"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: CHEBI:17515 phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    A phosphatidylglycerol has a glycerol backbone with two fatty acid esters and a phosphatidyl group (glycerol connected via phosphate).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of phosphorus (phosphate group)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphate group found"

    # Check for two ester groups (O-C=O)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2"

    # Check for glycerol backbone connected to phosphate
    # Main glycerol: two esters and a phosphate
    # Head glycerol: connected via phosphate, has hydroxyls
    # SMARTS pattern for glycerol-phosphate-glycerol structure
    pg_pattern = Chem.MolFromSmarts(
        "[CH2](-O-C(=O)*)-[CH](-OP(=O)([OX1])-[OX2]-[C@H]([CH2]O)[CH2]O)-[CH2](-O-C(=O)*)"
    )
    if not mol.HasSubstructMatch(pg_pattern):
        # Try without stereochemistry
        pg_pattern = Chem.MolFromSmarts(
            "[CH2](-O-C(=O)*)-[CH](-OP(=O)([OX1])-[OX2]-[CH]([CH2]O)[CH2]O)-[CH2](-O-C(=O)*)"
        )
        if not mol.HasSubstructMatch(pg_pattern):
            return False, "Glycerol-phosphate-glycerol backbone not found"

    # Check that the head glycerol has at least two hydroxyl groups
    head_glycerol_pattern = Chem.MolFromSmarts("[CH2]-[CH](-O)-[CH2]-O")
    if not mol.HasSubstructMatch(head_glycerol_pattern):
        return False, "Head glycerol missing hydroxyl groups"

    # Check fatty acid chains (at least 8 carbons each)
    # Approximate by checking molecular weight > 600 Da (varies, but PG is typically larger)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 600:
        return False, "Molecular weight too low for phosphatidylglycerol"

    return True, "Contains glycerol backbone with phosphatidyl group and two fatty acid esters"