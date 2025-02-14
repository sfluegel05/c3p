"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride has a glycerol backbone with a fatty acid chain at position 1 via an ester bond,
    and free hydroxyl groups at positions 2 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Glycerol backbone check (C-C-C with three oxygens and two hydroxyls and an ester)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]") # 
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"
    if len(glycerol_matches) > 1:
        return False, "More than 1 glycerol backbone found"
    
    glycerol_match = glycerol_matches[0] # get the tuple of atoms

    # Get the glycerol carbon atoms from the match, convert them to a list of integers
    glycerol_carbon_atoms = [int(atom) for atom in glycerol_match]

    #Check hydroxyl group at position 2
    glycerol_carbon_2_index = glycerol_carbon_atoms[1]
    glycerol_carbon_2 = mol.GetAtomWithIdx(glycerol_carbon_2_index)
    oxygen_neighbors_c2 = [neighbor for neighbor in glycerol_carbon_2.GetNeighbors() if neighbor.GetAtomicNum() == 8]
    if len(oxygen_neighbors_c2) != 1:
        return False, "Glycerol position 2 does not have one hydroxyl group"


    # Check that the first and third glycerol carbons have a hydroxyl group,
    # and for the first glycerol carbon we check it does NOT connect to more than 1 other carbon (meaning it will be connected to ester O)
    glycerol_carbon_1_index = glycerol_carbon_atoms[0]
    glycerol_carbon_1 = mol.GetAtomWithIdx(glycerol_carbon_1_index)
    oxygen_neighbors_c1 = [neighbor for neighbor in glycerol_carbon_1.GetNeighbors() if neighbor.GetAtomicNum() == 8]
    if len(oxygen_neighbors_c1) != 1:
      return False, "Glycerol position 1 does not have a hydroxyl group"
        
    glycerol_carbon_3_index = glycerol_carbon_atoms[2]
    glycerol_carbon_3 = mol.GetAtomWithIdx(glycerol_carbon_3_index)
    oxygen_neighbors_c3 = [neighbor for neighbor in glycerol_carbon_3.GetNeighbors() if neighbor.GetAtomicNum() == 8]
    if len(oxygen_neighbors_c3) != 1:
        return False, "Glycerol position 3 does not have a hydroxyl group"

    #2. Check for one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"
    
    ester_match = ester_matches[0]

    # Get the ester oxygen and carbonyl carbon atoms
    ester_oxygen_index = ester_match[0]
    ester_carbon_index = ester_match[1]
    ester_oxygen_atom = mol.GetAtomWithIdx(ester_oxygen_index)
    ester_carbon_atom = mol.GetAtomWithIdx(ester_carbon_index)
    
    # 4. Check if the ester is attached to the correct carbon of glycerol backbone (1-position)
    #Find if the ester oxygen is connected to the carbon at position 1 of the glycerol
    is_connected_to_glycerol_1 = False
    for neighbor in ester_oxygen_atom.GetNeighbors():
      if neighbor.GetIdx() == glycerol_carbon_1_index:
        is_connected_to_glycerol_1 = True
        break

    if not is_connected_to_glycerol_1:
      return False, "Ester is not attached to the 1-position of the glycerol backbone"
    
    #Check that this ester is ONLY attached to the glycerol backbone and not to anything else (besides fatty acid)
    glycerol_c_indexes = [glycerol_carbon_1_index, glycerol_carbon_2_index, glycerol_carbon_3_index]
    ester_oxygen_neighbors = [n.GetIdx() for n in ester_oxygen_atom.GetNeighbors()]
    
    ester_oxygen_neighbors_not_glycerol = [idx for idx in ester_oxygen_neighbors if idx not in glycerol_c_indexes ]
    if len(ester_oxygen_neighbors_not_glycerol) != 0:
        return False, "Ester is connected to something other than the glycerol backbone"

    
    #5. Check for fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") 
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "Missing fatty acid chain"
    
    # Ensure the fatty acid is connected to the carbonyl carbon of the ester
    is_connected_to_fatty_acid = False
    for neighbor in ester_carbon_atom.GetNeighbors():
      if neighbor.GetIdx() != ester_oxygen_index and neighbor.GetAtomicNum() == 6:
        is_connected_to_fatty_acid = True
        break

    if not is_connected_to_fatty_acid:
        return False, "Fatty acid chain not connected to the ester"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4: # Minimum of 4 rotatable bonds to check for reasonable fatty acid size
        return False, "Chain too short to be a fatty acid"

    # Check molecular weight - monoglycerides typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
      return False, "Molecular weight too low for monoglyceride"

    return True, "Contains a glycerol backbone with one fatty acid chain attached at the 1-position via an ester bond, and free hydroxyls at positions 2 and 3"