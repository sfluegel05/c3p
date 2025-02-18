"""
Classifies: CHEBI:16385 organic sulfide
"""
"""
Classifies: Organic sulfide (thioether)
Definition: Compounds having the structure R-S-R (with R ≠ H).
Such compounds were once called thioethers.
"""

from rdkit import Chem

def is_organic_sulfide(smiles: str):
    """
    Determines if a molecule is an organic sulfide (thioether) based on its SMILES string.
    An organic sulfide is defined as a molecule containing at least one R-S-R moiety,
    where both R groups are not hydrogen (typically organic substituents, e.g., carbon-based).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is an organic sulfide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a thioether group: carbon-sulfur-carbon.
    # Note: This pattern matches an S atom bonded to two carbon atoms.
    thioether_pattern = Chem.MolFromSmarts("[#6]-S-[#6]")
    
    # Check if the molecule has at least one substructure that matches the thioether pattern.
    if mol.HasSubstructMatch(thioether_pattern):
        return True, "Molecule contains at least one organic sulfide (thioether) group (R-S-R bond)"
    else:
        return False, "No organic sulfide (R-S-R) moiety found in the molecule"