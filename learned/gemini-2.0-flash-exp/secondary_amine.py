"""
Classifies: CHEBI:32863 secondary amine
"""
from rdkit import Chem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.
    A secondary amine is a nitrogen atom bonded to two carbon atoms and one hydrogen atom (implicitly or explicitly).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a secondary amine, where the nitrogen is bonded to 2 carbons.
    # The implicit hydrogen is accounted for because the nitrogen is explicitly defined to have 2 bonds and 1 implicit.
    secondary_amine_pattern = Chem.MolFromSmarts("[NHX2]([C])([C])")
    

    # Check for at least one match
    if not mol.HasSubstructMatch(secondary_amine_pattern):
        return False, "No secondary amine found (N with two carbons)"
    
    # Check if the nitrogen is quaternized
    quaternary_nitrogen_pattern = Chem.MolFromSmarts("[NX4]")
    if mol.HasSubstructMatch(quaternary_nitrogen_pattern):
         return False, "Nitrogen is quaternized, not a secondary amine"


    return True, "Molecule contains a secondary amine group"