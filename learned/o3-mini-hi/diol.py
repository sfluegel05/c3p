"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: diol (A compound that contains two hydroxy groups)

Note: The algorithm uses the pattern "[OX2H]" to find hydroxyl groups.
This pattern typically matches both alcoholic OH and the OH part of carboxylic acids.
However, since the definition of diol (as used here) is a compound that contains exactly two hydroxy groups,
we simply count the occurrences. Any molecule with exactly two -OH groups is classified as a diol.
For molecules with fewer or more than two hydroxyl groups the classification returns False.
"""

from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol (contains two hydroxy groups)
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a diol (exactly two hydroxy groups), False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a hydroxy group: oxygen atom with two connections and one hydrogen
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    
    # Find all matches of the hydroxy pattern in the molecule
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    num_hydroxy = len(hydroxy_matches)
    
    # Classify the molecule as diol if exactly two hydroxy groups are present
    if num_hydroxy == 2:
        return True, "Molecule contains exactly two hydroxy groups and is classified as a diol."
    else:
        return False, f"Molecule contains {num_hydroxy} hydroxy groups, which does not match the diol definition (exactly two required)."