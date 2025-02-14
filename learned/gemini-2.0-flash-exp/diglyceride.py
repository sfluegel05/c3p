"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: diglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    A diglyceride is a glycerol backbone with two fatty acid chains attached via ester bonds.
    The remaining hydroxyl group can be either free or alkylated

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diglyceride, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C)
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH][CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

     # Look for remaining hydroxyl or ether directly connected to glycerol backbone
    hydroxyl_or_ether_pattern = Chem.MolFromSmarts("[CH2X4,CHX4]-[OX2H0,OX2]-[#6]") #Check for -OH or -O-C, with -OH not having a hydrogen
    hydroxyl_or_ether_matches = mol.GetSubstructMatches(hydroxyl_or_ether_pattern)
    if len(hydroxyl_or_ether_matches) != 1:
            return False, f"Found {len(hydroxyl_or_ether_matches)} hydroxyl/ether groups attached to glycerol backbone, need 1"


    # Check for fatty acid chains attached to the two ester groups (long carbon chains)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])-[OX2]-[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"


    # Basic checks for number of carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 8:
        return False, "Too few carbons for diglyceride"
    if o_count < 4:
        return False, "Too few oxygens for diglyceride"
        

    return True, "Contains glycerol backbone with 2 fatty acid chains attached via ester bonds and a remaining hydroxyl or ether."