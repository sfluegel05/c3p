"""
Classifies: CHEBI:17297 UDP-sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    A UDP-sugar is a pyrimidine nucleotide-sugar having UDP as the nucleotide component
    attached to an unspecified sugar via an anomeric diphosphate linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): A tuple containing:
            - True if molecule is a UDP-sugar, False otherwise.
            - Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Identify the UDP core (relaxed stereochemistry)
    # Define a SMARTS pattern for the UDP core, without explicit Hs or charges.
    udp_core_pattern = Chem.MolFromSmarts('O=C1NC(=O)C=CN1C2C(O)C(O)C(COP(=O)(O)OP(=O)(O)O)C2O')
    udp_matches = mol.GetSubstructMatches(udp_core_pattern)
    if not udp_matches:
       return False, "No UDP core found"
    udp_match = udp_matches[0]

    # 2. Verify the diphosphate linkage
    # Check that the terminal phosphate of UDP is connected to a sugar moiety.
    # We'll check for a P-O-C bond where the C is not part of the UDP core
    linkage_pattern = Chem.MolFromSmarts("[#15](=[OX1])(-[OX2])-[CX4]")
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)

    if not linkage_matches:
       return False, "No diphosphate linkage found"

    #Check that the carbon is connected to the UDP via the oxygen.
    is_linked = False
    for match in linkage_matches:
        phosphorus_atom = mol.GetAtomWithIdx(match[0])
        oxygen_atom = mol.GetAtomWithIdx(match[1])
        carbon_atom = mol.GetAtomWithIdx(match[2])
        
        #Check if this carbon is not part of the ribose of UDP
        is_ribose_carbon = False
        for idx in udp_match:
            if idx == carbon_atom.GetIdx():
                is_ribose_carbon = True
        
        if not is_ribose_carbon:
            for neighbor in carbon_atom.GetNeighbors():
                if neighbor.GetIdx() == oxygen_atom.GetIdx():
                    is_linked=True
    
    if not is_linked:
        return False, "The diphosphate linkage is not connected to the sugar."

    # 3. Check the sugar moiety
    # Look for a carbon with an oxygen neighbor forming a ring system (or part of it)
    sugar_pattern = Chem.MolFromSmarts("[CX4][OX2]~[CX4]")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
      return False, "No sugar moiety found"

    
    #Check that this sugar is not part of the UDP ribose
    is_sugar_connected = False
    for sugar_match in sugar_matches:
        carbon_atom = mol.GetAtomWithIdx(sugar_match[0])
        is_ribose_carbon = False
        for idx in udp_match:
            if idx == carbon_atom.GetIdx():
                is_ribose_carbon = True
        if not is_ribose_carbon:
            is_sugar_connected = True
    
    if not is_sugar_connected:
      return False, "The sugar is part of the UDP core"

    # 4. Check the overall composition (relaxed constraints)
    # Check the correct number of atoms of each type
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)


    if n_count < 2:
      return False, "Too few nitrogens"
    if o_count < 8:
        return False, "Too few oxygens"
    if c_count < 8:
      return False, "Too few carbons"
    if p_count != 2:
      return False, "Must have exactly 2 phosphorus atoms"


    # If all checks passed, it's likely a UDP-sugar
    return True, "Contains UDP core with diphosphate linkage to a sugar"