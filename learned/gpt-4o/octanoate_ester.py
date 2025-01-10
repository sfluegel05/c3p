"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester contains the ester linkage where the carboxylic acid component is octanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Octanoic acid ester pattern: highlight ester linkage (C(=O)O) followed by 8-carbon chain
    # Use flexible branching that accounts for different positions and potential ring formations
    octanoate_ester_patterns = [
        Chem.MolFromSmarts("C(=O)OCCCCCCC[CH2]"),  # Linear octanoate ester
        Chem.MolFromSmarts("C(=O)O[C;R0]C(CCCCCC)"), # Branched ending or initial branching
        Chem.MolFromSmarts("[C;R]C(=O)OCCCCCCC"),    # Ring with ester
    ]

    # Check for the presence of the octanoic ester group among different patterns
    for pattern in octanoate_ester_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains octanoic acid ester linkage"
    
    return False, "Does not contain octanoic acid ester linkage"