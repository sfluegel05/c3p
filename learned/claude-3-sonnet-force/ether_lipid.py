"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: CHEBI:38114 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    Ether lipids contain at least one ether linkage (-O-) between a glycerol backbone
    and an alkyl chain. They may also contain ester linkages (-O-C(=O)-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for glycerol backbone pattern (C-C-C with 2-3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X2,CH2X3][CHX3,CHX4][CH2X2,CH2X3]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for ether linkages (-O-)
    ether_pattern = Chem.MolFromSmarts("[OX2]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if not ether_matches:
        return False, "No ether linkages found"
    
    # Check for long alkyl chains (>= 8 carbons)
    alkyl_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    alkyl_matches = mol.GetSubstructMatches(alkyl_pattern)
    if not alkyl_matches:
        return False, "No long alkyl chains found"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 12:
        return False, "Too few carbons for ether lipid"
    if o_count < 3:
        return False, "Too few oxygens for ether lipid"
    
    return True, "Contains glycerol backbone with at least one ether linkage and long alkyl chains"