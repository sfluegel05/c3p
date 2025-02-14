"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine is a six-membered ring with nitrogens at positions 1, 2, and 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2,4-triazine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the 1,2,4-triazine core with aromaticity enforcement
    triazine_core_smarts_aromatic = "n1-c=n-c=n-c1"
    triazine_pattern_aromatic = Chem.MolFromSmarts(triazine_core_smarts_aromatic)

    # Define SMARTS pattern for the 1,2,4-triazine core with possible non-aromatic form.
    triazine_core_smarts_nonaromatic = "N1-C-N-C=N-C1"
    triazine_pattern_nonaromatic = Chem.MolFromSmarts(triazine_core_smarts_nonaromatic)

    
    # Check if the molecule contains the core 1,2,4-triazine substructure (aromatic or non aromatic)
    if not mol.HasSubstructMatch(triazine_pattern_aromatic) and not mol.HasSubstructMatch(triazine_pattern_nonaromatic):
        return False, "Molecule does not contain the 1,2,4-triazine core."

    return True, "Molecule contains the 1,2,4-triazine core."