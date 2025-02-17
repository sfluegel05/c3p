"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: Organobromine compound
Definition: A compound containing at least one carbon–bromine bond.
This improved version uses a SMARTS matcher with a generic bond operator (~)
so that both aliphatic and aromatic C–Br bonds are captured.
"""

from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound is defined as a compound containing at least one carbon–bromine bond.
    
    This implementation uses a SMARTS pattern that matches any bond (using '~') between a 
    carbon atom (atomic number 6) and a bromine atom (atomic number 35). This allows capture
    of both aliphatic and aromatic bonds.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains at least one carbon–bromine bond, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern to match a bond (of any type) between a carbon (#6)
    # and bromine (Br). The operator '~' means "any bond" (single, aromatic, etc.)
    pattern = Chem.MolFromSmarts("[#6]~[Br]")
    if pattern is None:
        return False, "Error in SMARTS pattern"
    
    # Test if the molecule has any substructure that matches a carbon–bromine bond.
    if mol.HasSubstructMatch(pattern):
        return True, "Contains at least one carbon–bromine bond"
    else:
        return False, "No carbon–bromine bond found"

# Example usage (uncomment the following lines for testing):
# test_smiles = "Brc1ccccc1"  # bromobenzene
# result, reason = is_organobromine_compound(test_smiles)
# print(result, reason)