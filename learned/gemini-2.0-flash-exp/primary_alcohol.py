"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is a compound in which a hydroxy group, -OH, is attached to a
    saturated carbon atom which has either three hydrogen atoms attached to it or
    only one other carbon atom and two hydrogen atoms attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for primary alcohol carbons
    primary_alcohol_carbon_pattern1 = Chem.MolFromSmarts("[CH3][OH]")  # -CH2-OH
    primary_alcohol_carbon_pattern2 = Chem.MolFromSmarts("[CH2]([CX4])[OH]")  # -CH2-R-OH

    # Find matches for both patterns
    match1 = mol.GetSubstructMatches(primary_alcohol_carbon_pattern1)
    match2 = mol.GetSubstructMatches(primary_alcohol_carbon_pattern2)
    
    # Check if any of the patterns are present
    if not match1 and not match2:
      return False, "No primary alcohol group found"
      
    # Verify saturation of all carbon matches found
    for match in match1 + match2:
        # The carbon is the atom connected to oxygen, which is at index 0 or 1 depending on the pattern.
        if len(match) == 2:
           carbon_index = match[0]
        elif len(match) == 3:
           carbon_index = match[0]

        carbon_atom = mol.GetAtomWithIdx(carbon_index)
        if not carbon_atom.GetHybridization() == Chem.HybridizationType.SP3:
            return False, f"Carbon atom not saturated: atom {carbon_atom.GetIdx()}"

    return True, "Primary alcohol group found"