"""
Classifies: CHEBI:16385 organic sulfide
"""
from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide based on its SMILES string.
    An organic sulfide has the structure R-S-R', where R and R' are carbon-containing groups (not hydrogen).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for an organic sulfide
    sulfide_pattern = Chem.MolFromSmarts("[#6]~[#16]~[#6]")

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(sulfide_pattern):
        return True, "Has the structure R-S-R' (R and R' are carbon atoms)"
    else:
         # Check for presence of sulfur atom to make reason more precise
        sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'S']
        if not sulfur_atoms:
            return False, "No sulfur atom found"
        else:
            return False, "Sulfur is not bonded to at least two carbon atoms."