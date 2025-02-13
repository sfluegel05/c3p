"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.

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
    
    # Define patterns to search for diradylglycerol criteria
    # Glycerol backbone: A pattern with three carbons and attached oxygens could indicate a glycerol
    glycerol_pattern = Chem.MolFromSmarts("[C;X4H1]([O;X2])([C;X4H1])([C;X4H1])")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Ester linkage pattern: O=C-O
    ester_pattern = Chem.MolFromSmarts("O=C-O")
    ether_pattern = Chem.MolFromSmarts("C-O-C")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    
    # Ensure at least two substituent linkages (either ester or ether)
    if len(ester_matches) + len(ether_matches) < 2:
        return False, f"Expected at least 2 substituent groups, found {len(ester_matches) + len(ether_matches)}"
    
    # Count rotatable bonds in the substituent chains to ensure they are long enough
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Substituent chains too short to be considered valid diradylglycerol moieties"
    
    return True, "Contains glycerol backbone with two substituent groups attached via ester or ether bonds"