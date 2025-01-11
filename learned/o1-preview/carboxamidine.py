"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:37936 carboxamidine
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine is defined as 'Compounds having the structure RC(=NR)NR2.
    The term is used as a suffix in systematic nomenclature to denote the -C(=NH)NH2 group including its carbon atom.'

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the carboxamidine group
    carboxamidine_pattern = Chem.MolFromSmarts("[C;X3](=[N;X2])[N;X3]")

    # Search for the carboxamidine group in the molecule
    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains carboxamidine group (-C(=NR)NR2)"
    else:
        return False, "No carboxamidine group found"