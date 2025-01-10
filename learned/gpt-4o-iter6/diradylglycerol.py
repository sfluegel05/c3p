"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol has a glycerol backbone with two substituents connected 
    through ester or ether bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define glycerol backbone pattern: three carbons, oxygens at specific positions
    backbone_pattern = Chem.MolFromSmarts("[C;H2]([O;X2])[C;H]([O;X2])[C;H2][O;X2]")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No specific glycerol backbone found"
    
    # Define patterns for ester and ether linkages for chain identification
    ester_pattern = Chem.MolFromSmarts("O=C[O;!R1]")  # ester linkage, excluding ester inside rings
    ether_pattern = Chem.MolFromSmarts("[C]O[C;!R1]")  # ether linkage, excluding those in rings
    
    # Check for connections from the glycerol backbone
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)

    # Verify exactly two substituent attachments via ester or ether bonds from glycerol
    if len(ester_matches) + len(ether_matches) != 2:
        return False, f"Expected exactly 2 substituent groups, found {len(ester_matches) + len(ether_matches)}"

    # Assess length and nature of attaching chains to ensure valid substituent chains
    for match in ester_matches + ether_matches:
        carbon_chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic() is False)
        if carbon_chain_length < 7:  # arbitrary long-chain threshold
            return False, "Substituent chains too short to be considered valid diradylglycerol moieties"
    
    return True, "Contains a glycerol backbone with exactly two substituent groups connected via ester or ether bonds"