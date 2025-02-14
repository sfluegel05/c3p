"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid has a glycerol backbone, a phosphate group, and two fatty acid chains.
    One fatty acid can be linked to the phosphate.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone with phosphate attached
    # Here we enforce that one oxygen of the phosphate is directly attached to the glycerol
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4]([OX2][P](=[OX1])([OX2])([OX2]))[CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_phosphate_pattern):
        return False, "No glycerol backbone with attached phosphate found"

    # Look for 2 ester groups attached to the glycerol backbone or directly to phosphate group
    # The ester must be linked to the glycerol backbone
    ester_pattern1 = Chem.MolFromSmarts("[CH2X4,CHX4]-O-[CX3](=[OX1])") # ester to glycerol
    ester_pattern2 = Chem.MolFromSmarts("[OX2][P](=[OX1])([OX2])([OX2])-[OX2][CX3](=[OX1])") # ester to phosphate
    
    ester_matches1 = mol.GetSubstructMatches(ester_pattern1)
    ester_matches2 = mol.GetSubstructMatches(ester_pattern2)
    
    if len(ester_matches1) + len(ester_matches2) < 2 :
        return False, f"Found only {len(ester_matches1) + len(ester_matches2)} ester groups attached to the glycerol or phosphate, need at least 2"


    #Check for long fatty acid chains attached to glycerol via ester or directly attached to phosphate via ester
    fatty_acid_pattern_ester = Chem.MolFromSmarts("[CH2X4,CHX4]-O-[CX3](=[OX1])-[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_pattern_phosphate = Chem.MolFromSmarts("[OX2][P](=[OX1])([OX2])([OX2])-[OX2][CX3](=[OX1])-[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    
    fatty_acid_matches_ester = mol.GetSubstructMatches(fatty_acid_pattern_ester)
    fatty_acid_matches_phosphate = mol.GetSubstructMatches(fatty_acid_pattern_phosphate)

    if len(fatty_acid_matches_ester) + len(fatty_acid_matches_phosphate) < 2:
        return False, f"Missing fatty acid chains, found {len(fatty_acid_matches_ester) + len(fatty_acid_matches_phosphate)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"

    # Count carbons and oxygens - expect at least 20 carbons and at least 6 oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 20:
        return False, "Too few carbons for a phosphatidic acid"
    if o_count < 6:
        return False, "Too few oxygens for a phosphatidic acid"


    return True, "Contains glycerol backbone with a phosphate group and two fatty acid chains attached via ester bonds"