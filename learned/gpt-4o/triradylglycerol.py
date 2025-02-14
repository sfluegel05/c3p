"""
Classifies: CHEBI:76579 triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify the glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[C@H](O)[C@H](O)[C@H](O)")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Identify the ester bonds to long chains [-C(=O)O- or -CO-]
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C@H]")
    acyl_pattern = Chem.MolFromSmarts("C(=O)O")
    alkyl_pattern = Chem.MolFromSmarts("CO")
    alk1enyl_pattern = Chem.MolFromSmarts("CO[C@H]=C")

    matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check for ester linkages in all 3 positions
    if len(matches) < 3:
        return False, f"Found only {len(matches)} ester/alkyl/alk-1-enyl linkages, need all 3"

    # Verify diversity of substituent types
    acyl_match_count = len(mol.GetSubstructMatches(acyl_pattern))
    alkyl_match_count = len(mol.GetSubstructMatches(alkyl_pattern))
    alk1enyl_match_count = len(mol.GetSubstructMatches(alk1enyl_pattern))

    if acyl_match_count + alkyl_match_count + alk1enyl_match_count < 3:
        return False, "Less than 3 valid acyl/alkyl/alk-1-enyl substituents found"

    return True, "Is a triradylglycerol with appropriate substituents at all three positions"