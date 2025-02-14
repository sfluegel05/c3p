"""
Classifies: CHEBI:37739 glycerophospholipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is a glycerol backbone with two fatty acid chains and a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 3 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
       return False, "No glycerol backbone found"


    # Look for phosphate group connected to glycerol through an ester bond
    phosphate_pattern = Chem.MolFromSmarts("[CX4][OX2][PX4](=[OX1])([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
      return False, "No phosphate group attached to glycerol via an ester bond found"


    # Look for two ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
      return False, "Less than two ester groups found"


    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Missing fatty acid chains, found {len(fatty_acid_matches)}"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 6:
        return False, "Chains too short to be fatty acids"

     # Check the connectivity to ensure phosphate is connected to glycerol terminal carbon
    glycerol_carbons = mol.GetSubstructMatches(glycerol_pattern)
    
    if glycerol_carbons and phosphate_matches:
        glycerol_carbon_atoms = [mol.GetAtomWithIdx(i) for i in glycerol_carbons[0]]
        phosphorous_atoms = [mol.GetAtomWithIdx(i[2]) for i in phosphate_matches]
        
        if not any (gc.GetIdx() == pa.GetNeighbors()[0].GetIdx() for gc in [glycerol_carbon_atoms[0], glycerol_carbon_atoms[2]] for pa in phosphorous_atoms):
            return False, "Phosphate not attached to a terminal glycerol carbon"
        
    return True, "Contains glycerol backbone with 2 fatty acid chains and a phosphate group linked via an ester bond."