"""
Classifies: CHEBI:18035 diglyceride
"""
from rdkit import Chem

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    A diglyceride is a glycerol molecule with two of its hydroxy groups acylated.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a diglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify the glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) == 0:
        return False, "No glycerol backbone found"
    
    # Identify the ester linkage pattern (RCO-OC)
    ester_pattern = Chem.MolFromSmarts("C(=O)O[CH2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check that there are exactly two ester linkages
    if len(ester_matches) != 2:
        return False, f"Found {len(ester_matches)} ester linkages, require exactly 2 for diglyceride"
    
    # Check if each ester linkage is connected to the glycerol backbone
    glycerol_carbon_indices = set()
    for match in glycerol_matches:
        glycerol_carbon_indices.update(match[1:])  # The indices of the carbons in glycerol
    
    ester_connected_to_glycerol = 0
    for ester_match in ester_matches:
        ester_o = ester_match[-1]  # The last carbon is the one connecting to glycerol
        ester_carbon_candidate = mol.GetAtomWithIdx(ester_o).GetNeighbors()
        if any(neighbor.GetIdx() in glycerol_carbon_indices for neighbor in ester_carbon_candidate):
            ester_connected_to_glycerol += 1
    
    if ester_connected_to_glycerol != 2:
        return False, "Ester groups not properly linked to two positions on the glycerol backbone"
    
    return True, "Contains glycerol backbone with two fatty acid chains attached via ester bonds"