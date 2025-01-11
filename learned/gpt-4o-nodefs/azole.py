"""
Classifies: CHEBI:68452 azole
"""
"""
Classifies: azole
"""
from rdkit import Chem

def is_azole(smiles: str):
    """
    Determines if a molecule contains an azole ring based on its SMILES string.
    Azoles are characterized by a 5-membered heterocyclic ring containing nitrogen atoms and
    optionally other heteroatoms like oxygen or sulfur.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an azole, False otherwise
        str: Reason for classification
    """

    # Revised SMARTS patterns for different types of azoles
    azole_patterns = {
        "pyrazole": "[nH]1nccc1",
        "imidazole": "n1c[nH]cc1",
        "oxazole": "o1cncc1",
        "thiazole": "c1cscn1",
        "1,2,3-triazole": "n1ncnn1",
        "1,2,4-triazole": "n1ncnc1",
        "isoxazole": "o1ncc[c,n]1",
        "benzotriazole": "c1ncn[nH]c2ccccc12"
    }

    # Parse SMILES string into molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of any azole substructure
    for name, pattern in azole_patterns.items():
        azole_substruct = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(azole_substruct):
            return True, f"Molecule contains {name} ring"

    return False, "No azole ring found"