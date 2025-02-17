"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: alkanesulfonate oxoanion
Definition: An alkanesulfonate in which the carbon at position 1 is attached to R, 
which can represent hydrogens, a carbon chain, or other groups.
This program checks for the presence of an sp3 (non‐aromatic) carbon directly attached to 
a sulfonate group S(=O)(=O)[O-] in the molecule.
"""

from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion is defined here as having a –S(=O)(=O)[O-] group directly attached
    to an sp3 (non-aromatic) carbon.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains an alkanesulfonate oxoanion moiety, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for an alkanesulfonate group.
    # The pattern "[C;!a]S(=O)(=O)[O-]" matches:
    # • A carbon atom (C) that is not aromatic (;!a),
    # • Directly bonded to a sulfur atom (S) with two double-bonded oxygens (=O) 
    #   and one anionic oxygen ([O-]).
    sulfonate_pattern = Chem.MolFromSmarts("[C;!a]S(=O)(=O)[O-]")
    
    # Search for the pattern in the molecule
    matches = mol.GetSubstructMatches(sulfonate_pattern)
    if not matches:
        return False, "No alkanesulfonate group (C–S(=O)(=O)[O-]) found"
    
    # If at least one match is found we consider the molecule to have
    # an alkanesulfonate oxoanion moiety.
    return True, "Contains an alkanesulfonate oxoanion moiety"