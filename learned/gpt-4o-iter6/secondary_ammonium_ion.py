"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
from rdkit import Chem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is an organic cation obtained by protonation of a secondary amine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
   
    # Parse the SMILES string to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced SMARTS pattern for secondary ammonium ion
    # Nitrogen positively charged [NH2+] bound to exactly two carbon atoms [C], no additional neighbors
    pattern = Chem.MolFromSmarts("[NH2+][C][C]")

    # Check if the molecule matches the secondary ammonium ion pattern
    matches = mol.GetSubstructMatches(pattern)
    
    if matches:
        return True, "Matches the secondary ammonium ion structure pattern"
    else:
        return False, "Does not match the secondary ammonium ion structure pattern"