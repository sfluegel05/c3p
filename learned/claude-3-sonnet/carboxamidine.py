"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:35813 carboxamidine

Carboxamidines are compounds having the structure RC(=NR)NR2. The term is used
as a suffix in systematic nomenclature to denote the -C(=NH)NH2 group including
its carbon atom.
"""

from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxamidine substructure
    carboxamidine_pattern = Chem.MolFromSmarts("[NX3][CX3]([NX3])=[NX2]")
    matches = mol.GetSubstructMatches(carboxamidine_pattern)

    if matches:
        return True, "Molecule contains the carboxamidine substructure R-C(=NR)-NR2"
    else:
        return False, "No carboxamidine substructure found"