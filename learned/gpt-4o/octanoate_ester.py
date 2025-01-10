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
    
    # Define detailed octanoic acid ester pattern: ester group (O=C-O) followed by a straight or branched 8-carbon chain
    octanoate_ester_pattern = Chem.MolFromSmarts("C(=O)OCCCCCCCC")
    branched_octanoate_pattern = Chem.MolFromSmarts("C(=O)O[C;R0]1CCCCCC1")  # Anchoring for ring closure or branching at first carbon
    
    # Check for the presence of the octanoic ester group with branching consideration
    if mol.HasSubstructMatch(octanoate_ester_pattern) or mol.HasSubstructMatch(branched_octanoate_pattern):
        return True, "Contains octanoic acid ester linkage"
    
    return False, "Does not contain octanoic acid ester linkage"