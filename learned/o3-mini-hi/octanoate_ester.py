"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: Octanoate Ester
Any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid),
i.e. the acyl portion is CH3(CH2)6C(=O)O–.
"""

from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    
    An octanoate ester is defined as any fatty acid ester where the acyl group 
    originates from octanoic acid. This means the acyl group should correspond to:
    CH3(CH2)6C(=O)O–
    
    The function searches for an exact substructure pattern that matches:
         CH3 - CH2 - CH2 - CH2 - CH2 - CH2 - CH2 - C(=O)O
    which guarantees that the acyl chain has the proper number of carbons and that the
    terminal methyl group and the carbonyl system are in the correct positions.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a proper octanoate ester substructure is present, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the explicit SMARTS pattern for an octanoate ester acyl group:
    # CH3-CH2-CH2-CH2-CH2-CH2-CH2-C(=O)O
    octanoate_smarts = "[CH3]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-[CH2]-C(=O)O"
    octanoate_pattern = Chem.MolFromSmarts(octanoate_smarts)
    if octanoate_pattern is None:
        return False, "Internal error: invalid SMARTS pattern"
    
    # Check if the molecule has the octanoate ester substructure.
    if mol.HasSubstructMatch(octanoate_pattern):
        return True, "Contains an octanoate ester group (acyl derived from octanoic acid)."
    
    return False, "No octanoate ester group found."