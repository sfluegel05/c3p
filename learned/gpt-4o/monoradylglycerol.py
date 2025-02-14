"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol has a glycerol backbone with one acyl, alkyl, or alk-1-enyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")

    # Check for exactly one glycerol backbone match
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) != 1:
        return False, "Glycerol backbone not found or multiple present"

    # Define patterns to recognize substituents
    acyl_pattern = Chem.MolFromSmarts("C(=O)[O;R]")  # R keeps it connected as part of larger structure
    alkyl_pattern = Chem.MolFromSmarts("[C;R]([C;R])[C;R]([C;R])[C;R]")  # at least 4 connected carbons
    alk1enyl_pattern = Chem.MolFromSmarts("C=C")  # ensure presence of double bonds

    # Check for substituents
    acyl_match = mol.HasSubstructMatch(acyl_pattern)
    alkyl_match = mol.HasSubstructMatch(alkyl_pattern)
    alk1enyl_match = mol.HasSubstructMatch(alk1enyl_pattern)

    # Count valid substituents
    num_substituents = sum([acyl_match, alkyl_match, alk1enyl_match])

    # Ensure there is exactly one substituent
    if num_substituents == 1:
        return True, "Contains a glycerol backbone with a single acyl, alkyl or alk-1-enyl substituent"
    
    return False, "Glycerol backbone does not have exactly one suitable substituent"