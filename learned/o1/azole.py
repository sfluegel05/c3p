"""
Classifies: CHEBI:68452 azole
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    An azole is defined as any monocyclic heteroarene consisting of a five-membered ring containing nitrogen.
    Azoles can also be part of fused ring systems and may contain other heteroatoms like nitrogen, sulfur, or oxygen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an azole, False otherwise
        str: Reason for classification
    """

    # Parse SMILES and sanitize molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Update aromaticity information
    Chem.SanitizeMol(mol)
    Chem.SetAromaticity(mol, Chem.AromaticityModel.AROMATICITY_MMFF94)

    # Define SMARTS patterns for azoles
    azole_smarts = [
        'n1cccn1',  # Pyrazole
        'n1cncc1',  # Imidazole
        'n1cocc1',  # Oxazole
        'n1csc(c1)',  # Thiazole
        'n1c[nH]nc1',  # Triazole
        'n1nnnn1',  # Tetrazole
        '[nH]1cccc1',  # Pyrrole
        'n1c[nH]cc1',  # Indazole
        'n1cccc1',  # Generic azole
        'n1ccnc1',  # Isoxazole
        'n1csc1',  # Isothiazole
    ]

    found_azole = False
    reasons = []

    for pattern in azole_smarts:
        azole_mol = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(azole_mol):
            found_azole = True
            # Get matching atoms to provide a reason
            matches = mol.GetSubstructMatches(azole_mol)
            reasons.append(f"Contains azole ring matched by SMARTS pattern '{pattern}'")
            break

    if found_azole:
        return True, reasons[0]
    else:
        return False, "No azole ring found (not an azole)"