"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine is defined as having the structure RC(=NR)NR2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # This SMARTS pattern accounts for the carboxamidine group: RC(=NR)NR2
    # It checks for the central carbon with two nitrogens attached
    # and allows for varied R groups
    carboxamidine_pattern = Chem.MolFromSmarts('C(=[N;H1,H2,H3])(N)[N;!H0]')
    
    # Check for matches
    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains carboxamidine (RC(=NR)NR2) group"
    
    return False, "No carboxamidine structure found"

# Example Usage:
print(is_carboxamidine("CC(Oc1c(Cl)cccc1Cl)C1=NCCN1"))  # Expected: True for lofexidine