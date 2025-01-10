"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine encompasses structures with a configuration RC(=NR)NR2,
    accounting for various R groups and structural contexts.

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
    
    # Refined SMARTS pattern for carboxamidine
    # This pattern requires a carbon double-bonded to a nitrogen, which is bonded
    # to two additional substituents which can vary (e.g., hydrogen, carbon chains)
    # Adjusting the SMARTS to better capture typical carboxamidine structures:
    carboxamidine_pattern = Chem.MolFromSmarts('[C;H0;$(C(=N)[N&R]);$(C-N)][NX3][R]')
    
    # Check for matches
    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains carboxamidine (RC(=NR)NR2) group"
    
    return False, "No carboxamidine structure found"