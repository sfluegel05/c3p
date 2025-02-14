"""
Classifies: CHEBI:76575 monoradylglycerol
"""
from rdkit import Chem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol has a glycerol backbone bearing a single acyl, alkyl, or alk-1-enyl substituent.

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

    # Define patterns for a glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts("C(CO)CO")
    
    # Check for the presence of a glycerol backbone
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Count ester/hydroxy linkages (acyl group check)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Count the number of long carbon chains (assuming a carbon chain of length >= 4 is a substituent)
    long_chain_pattern = Chem.MolFromSmarts("C(C(C(C)C)C)C")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)

    # Check for the presence of a single acyl/alkyl/alk-1-enyl substituent
    # Note: For strict monoradylglycerol checking, we assume one ester bond and no extra long chains
    if len(ester_matches) == 1 or len(long_chain_matches) == 1:
        return True, "Contains a glycerol backbone with a single acyl, alkyl or alk-1-enyl substituent"
    
    # If none of the patterns match specific criteria, it is not a monoradylglycerol
    return False, "Glycerol backbone does not have exactly one suitable substituent"