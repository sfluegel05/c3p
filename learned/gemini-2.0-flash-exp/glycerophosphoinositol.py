"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol has a glycerol backbone with two fatty acid chains and a phosphate group,
    where the phosphate is esterified with inositol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # 2. Check for the presence of a phosphate group attached to the glycerol
    phosphate_glycerol_pattern = Chem.MolFromSmarts("[CH2X4][OX2][P](=[OX1])([OX2])[OX2]")
    phosphate_glycerol_pattern2 = Chem.MolFromSmarts("[CHX4][OX2][P](=[OX1])([OX2])[OX2]")
    
    if not (mol.HasSubstructMatch(phosphate_glycerol_pattern) or mol.HasSubstructMatch(phosphate_glycerol_pattern2)):
          return False, "No phosphate group linked to glycerol found"
    
    # 3. Check for the presence of inositol attached to phosphate
    inositol_phosphate_pattern = Chem.MolFromSmarts("[P]([OX2])([OX2])-[OX2]-[C]1([O])[C]([O])[C]([O])[C]([O])[C]([O])[C]([O])1")
    if not mol.HasSubstructMatch(inositol_phosphate_pattern):
        return False, "No inositol group linked to the phosphate group found"

    # 4. Check for exactly two ester groups
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # 5. Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"

    # 6. Check for exactly one phosphate group attached to the glycerol
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
       return False, f"Found {p_count} phosphate groups, need exactly 1"

    return True, "Contains glycerol backbone, 2 fatty acid chains, inositol, and a phosphate group attached to inositol and glycerol"