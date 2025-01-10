"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxamidine pattern C(=NR)NR2
    # Create a query for C(=N)N structure
    carboxamidine_pattern = Chem.MolFromSmarts("C(=N)N")
    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains carboxamidine (C(=NR)NR2) group"

    # If no matches found, it's not a carboxamidine
    return False, "No carboxamidine structure found"

# Example Usage
print(is_carboxamidine("CC(Oc1c(Cl)cccc1Cl)C1=NCCN1"))  # Should be True for lofexidine