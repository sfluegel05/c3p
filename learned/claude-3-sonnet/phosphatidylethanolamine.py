"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: CHEBI:17355 phosphatidylethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine (PE) based on its SMILES string.
    A PE is a glycerophospholipid with a phosphatidyl group esterified to the hydroxy group of ethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a PE, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphate group (-O-P(=O)(-O)(-O)-) attached to glycerol
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=[OX1])([OX2])[OX2]")
    phosphate_match = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_match:
        return False, "No phosphate group found"
    
    # Look for ethanolamine group (-O-CH2-CH2-NH2) attached to phosphate
    ethanolamine_pattern = Chem.MolFromSmarts("[OX2]CCN")
    ethanolamine_match = mol.GetSubstructMatches(ethanolamine_pattern)
    if not ethanolamine_match:
        return False, "No ethanolamine group found"
    
    # Look for 2 ester groups (-O-C(=O)-) attached to glycerol
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"
    
    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, "Too few carbons for phosphatidylethanolamine"
    if o_count != 7:
        return False, "Must have exactly 7 oxygens (2 esters, 1 phosphate, 1 ethanolamine)"
    
    return True, "Contains glycerol backbone with 2 fatty acid chains, phosphate group, and ethanolamine group"