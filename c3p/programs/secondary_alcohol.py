"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: Secondary Alcohol
Definition:
    A secondary alcohol is a compound in which a hydroxy group (-OH) is attached to a
    saturated carbon atom which has two other carbon atoms attached to it.
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    
    A secondary alcohol has an -OH group attached to an sp3 carbon that is bonded to two 
    carbon atoms (and one hydrogen, given that the carbon is tetrahedral).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a secondary alcohol:
    # "[C;X4;H1]([#6])([#6])[OH]"
    # Explanation:
    #   [C;X4;H1] : a tetrahedral (sp3) carbon with exactly one hydrogen (implying a secondary center)
    #   ([#6])([#6]) : the carbon must be bonded to two carbon atoms
    #   [OH] : that same carbon is bonded to an -OH group.
    sec_alcohol_smarts = "[C;X4;H1]([#6])([#6])[OH]"
    pattern = Chem.MolFromSmarts(sec_alcohol_smarts)
    if pattern is None:
        return False, "Error processing SMARTS pattern for secondary alcohol"
    
    # Search for the secondary alcohol substructure in the molecule
    if not mol.HasSubstructMatch(pattern):
        return False, "No secondary alcohol functional group found"
    
    # Optionally, you could retrieve and inspect matches:
    matches = mol.GetSubstructMatches(pattern)
    if len(matches) == 0:
        return False, "No matching secondary alcohol groups found"
    
    return True, "Contains at least one secondary alcohol group"