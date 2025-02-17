"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: Octanoate ester (any fatty acid ester in which the carboxylic acid component is octanoic acid)
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is defined as any fatty acid ester where the acid component 
    is octanoic (caprylic) acid, i.e. contains an octanoate (CH3(CH2)6C(=O)O-) group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if an octanoate ester functional group is present, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern for an octanoate ester group.
    # Octanoic acid is CH3(CH2)6COOH, and its ester has the acyl component: CH3(CH2)6C(=O)O-
    # In SMARTS, this can be represented as: "CCCCCCCC(=O)O"
    # Explanation: "CCCCCCCC" represents a chain of 8 carbons where the first 7 are the alkyl chain
    # and the 8th a carbonyl carbon (C(=O)O part) attached to the ester oxygen.
    octanoate_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)O")
    if octanoate_pattern is None:
        return False, "Internal error: invalid SMARTS pattern"
    
    # Check if the molecule contains the octanoate ester substructure.
    if mol.HasSubstructMatch(octanoate_pattern):
        return True, "Contains an octanoate ester group (acyl derived from octanoic acid)"
    else:
        return False, "No octanoate ester group found"