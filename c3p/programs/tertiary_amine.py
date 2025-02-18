"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: tertiary amine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is a nitrogen atom bonded to three carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create SMARTS for tertiary amine (N connected to 3 carbons)
    tertiary_amine_pattern = Chem.MolFromSmarts("[NX3]([CX4])([CX4])[CX4]")

    # Check for substructure match
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Tertiary amine"
    else:
        return False, "Not a tertiary amine"