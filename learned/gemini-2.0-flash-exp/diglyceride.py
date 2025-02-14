"""
Classifies: CHEBI:18035 diglyceride
"""
"""
Classifies: diglyceride
"""
from rdkit import Chem
from rdkit.Chem import AllChem
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

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached via single bonds)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4](-[OX2])-[CHX4](-[OX2])-[CH2X4](-[OX2])")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for 2 ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 2:
      return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    # Check for a remaining alcohol or ether on glycerol backbone:
    hydroxyl_or_ether_pattern = Chem.MolFromSmarts("[CH2X4,CHX4](-[OX2]) | [CH2X4,CHX4,CHX3]-[OX2]-[CX4]")
    hydroxyl_or_ether_matches = mol.GetSubstructMatches(hydroxyl_or_ether_pattern)

    glycerol_carbons_pattern = Chem.MolFromSmarts("[CH2X4](-[OX2])-[CHX4](-[OX2])-[CH2X4](-[OX2])")
    glycerol_carbon_matches = mol.GetSubstructMatches(glycerol_carbons_pattern)

    if len(glycerol_carbon_matches) != 1:
        return False, "Incorrect number of glycerol backbones"

    if len(hydroxyl_or_ether_matches) < 1:
        return False, "No remaining hydroxyl or ether on glycerol backbone"
    
    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3](=[OX1])-[OX2]-[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, got {len(fatty_acid_matches)}"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - diglycerides typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for diglyceride"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 10:
        return False, "Too few carbons for diglyceride"

    if o_count < 5:
        return False, "Too few oxygens for diglyceride"

    return True, "Contains glycerol backbone with 2 fatty acid chains attached via ester bonds"