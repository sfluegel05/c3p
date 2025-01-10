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

    # Enhanced SMARTS patterns for various azoles
    # These patterns aim to capture standard and tautomeric forms
    azole_patterns = {
        "pyrazole": "[nH]1ncc[nH]1",
        "imidazole": "c1cn[nH]c1",
        "oxazole": "o1cncn1",
        "thiazole": "c1cscn1",
        "1,2,3-triazole": "n1nncn1",
        "1,2,4-triazole": "n1ncnc1",
        "isoxazole": "[o]1[nH]c[nH]c1",
        "benzotriazole": "c1nncc2c1cccc2",
        "pyrazolone": "C1=C(C=O)N[H]N1",
        "fused_azole": "c1[nH]nnc2ccccc12"  # Example of a fused system
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