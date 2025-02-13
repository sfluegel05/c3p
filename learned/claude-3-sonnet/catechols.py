"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: CHEBI:27027 catechol
Any compound containing an o-diphenol component
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catechol(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol is any compound containing an o-diphenol component (two hydroxyl groups
    attached to adjacent carbon atoms, which may or may not be part of a ring system).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enumerate tautomers
    tautomers = [mol]
    try:
        tautomers = list(AllChem.EnumerateCanonicalTautomers(mol))
    except ValueError:
        pass

    # Check for adjacent hydroxyl groups
    for tautomer in tautomers:
        for bond in tautomer.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetSymbol() == 'O' and atom1.GetIsAromatic() and atom2.GetSymbol() == 'C':
                for neighbor in atom2.GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and neighbor.GetIsAromatic():
                        return True, "Contains an o-diphenol component"

    return False, "Does not contain an o-diphenol component"