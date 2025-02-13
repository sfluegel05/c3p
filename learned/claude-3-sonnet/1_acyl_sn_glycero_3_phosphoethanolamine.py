"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: CHEBI:17901 1-acyl-sn-glycero-3-phosphoethanolamine
A 1-O-acylglycerophosphoethanolamine having (R)-configuration.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Look for phosphoethanolamine group (-O-P(=O)(-O-)-O-CH2-CH2-NH2)
    pe_pattern = Chem.MolFromSmarts("OP(OCCN)(=O)O[CX4H2][CX4H2]O")
    if not mol.HasSubstructMatch(pe_pattern):
        return False, "No phosphoethanolamine group found"

    # Check for (R)-configuration at the chiral carbon
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnspec=True)
    if len(chiral_centers) != 1:
        return False, "Expected exactly one chiral center"

    chiral_atom = mol.GetAtomWithIdx(chiral_centers[0][0])
    if chiral_atom.GetProp('_CIPCode') != 'R':
        return False, "Chiral center not in (R)-configuration"

    # Additional checks for long fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, f"Missing fatty acid chain, got {len(fatty_acid_matches)}"

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Fatty acid chain too short"

    return True, "Molecule matches the structure of a 1-acyl-sn-glycero-3-phosphoethanolamine"