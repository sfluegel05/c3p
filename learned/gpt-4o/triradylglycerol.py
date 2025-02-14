"""
Classifies: CHEBI:76579 triradylglycerol
"""
from rdkit import Chem

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol is characterized by a glycerol backbone with acyl, alkyl, or alk-1-enyl substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, "Invalid SMILES string"

        # Pattern to identify glycerol backbone with three possible substituents
        glycerol_pattern = Chem.MolFromSmarts("OCC(O)C(O)")  # Glycerol pattern
        matching_atoms = mol.GetSubstructMatches(glycerol_pattern)
        if not matching_atoms:
            return False, "No glycerol backbone found"

        # Identify ester or ether linkages in the molecule
        ester_pattern = Chem.MolFromSmarts("C(=O)O[CX4]")
        ether_pattern = Chem.MolFromSmarts("CO[CX4]")

        ester_matches = mol.GetSubstructMatches(ester_pattern)
        ether_matches = mol.GetSubstructMatches(ether_pattern)
        
        # Triradylglycerol should have three linkages combining esters and ethers
        if len(ester_matches) + len(ether_matches) != 3:
            return False, "Does not have three ester or ether groups attached to glycerol backbone"

        # Check for substituent diversity based on acyl, alkyl, alk-1-enyl
        substituent_count = len(set(ester_matches + ether_matches))
        if substituent_count < 3:
            return False, "Insufficient substituent diversity"

        return True, "Is a triradylglycerol with appropriate substituents"

    except Exception as e:
        return None, f"Exception occurred: {str(e)}"