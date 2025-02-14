"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin is a phosphatidylglycerol composed of two molecules of phosphatidic acid
    covalently linked to a molecule of glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Central glycerol backbone
    central_glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(central_glycerol_pattern):
        return False, "No central glycerol backbone found"
        
    # Check for glycerol connected to a phosphate group
    glycerol_phosphate_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4][OX2][P](=[OX1])([OX2])([OX2])")
    glycerol_phosphate_matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
    if len(glycerol_phosphate_matches) < 2:
      return False, f"Found {len(glycerol_phosphate_matches)} glycerol-phosphate groups, need at least 2"
      
    #Check for fatty acid chains attached to the glycerol in the phosphatidic group
    fatty_acid_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4][OX2][CX3](=[OX1])~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 4:
       return False, f"Found {len(fatty_acid_matches)} fatty acid chains connected to glycerol, need at least 4"

    # Check that there are 4 ester groups
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 4:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 4"


    # Check for the number of glycerol carbons: central glycerol + two from each phosphatidic acid.
    # Each phosphatidic acid should have at least 3 carbons from glycerol, so total glycerol carbons should be at least 3 + (2 *3 ) = 9
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 9: 
        return False, f"Too few carbons: {c_count}, less than expected"
    
    # Check number of phosphorous atoms (should be 2)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 2:
        return False, f"Must have exactly 2 phosphorus atoms, found {p_count}"

    # Check for the number of oxygens (should be more than 12)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 12:
        return False, f"Too few oxygens: {o_count}, less than expected"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 15:
        return False, "Chains too short to be fatty acids"

    # Check molecular weight - cardiolipins typically >1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight too low for cardiolipin"

    return True, "Contains two phosphatidic acids linked to a glycerol backbone"