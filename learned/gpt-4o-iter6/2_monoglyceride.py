"""
Classifies: CHEBI:17389 2-monoglyceride
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride is defined by an ester linkage at the second carbon of a glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for 2-monoglyceride pattern with ester at the second carbon
    glycerol_2_mono_pattern = Chem.MolFromSmarts("C(CO)C(OC(=O))CO")
    
    # Validate that the ester is specifically linked to the second carbon
    matches = mol.GetSubstructMatches(glycerol_2_mono_pattern)
    if not matches:
        return False, "No specific 2-monoglyceride pattern found with esterification at C2"

    # Check the matches for the correct configuration
    for match in matches:
        c1_index, c2_index, ester_oxygen_index, c3_index = match
        
        # Ensure C1 and C3 have the hydroxyl group ([CX2H] corresponds to a carbon with one hydrogen and -OH)
        c1_oh_check = mol.GetAtomWithIdx(c1_index).GetDegree() == 3  # Non-terminal carbon, typically connected to -OH 
        c3_oh_check = mol.GetAtomWithIdx(c3_index).GetDegree() == 3  # Same for C3

        # Ester oxygen should be bonded to C2
        ester_oxygen_bonded_to_c2 = any(bond.GetEndAtomIdx() == c2_index for bond in mol.GetAtomWithIdx(ester_oxygen_index).GetBonds())

        if c1_oh_check and c3_oh_check and ester_oxygen_bonded_to_c2:
            return True, "Contains 2-monoglyceride structure with acyl group esterified at the second position"

    return False, "The structure does not match the specific pattern of a 2-monoglyceride"