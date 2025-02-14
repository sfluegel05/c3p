"""
Classifies: CHEBI:28802 flavonols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a flavone with a hydroxyl group at the 3-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavonol core pattern using SMARTS. It specifies:
    # - Aromatic ring with substitutions, connected to
    # - The heterocyclic ring with oxygen and a carbonyl group
    # - A hydroxyl group directly connected to position 3 of the heterocyclic ring
    flavonol_pattern = Chem.MolFromSmarts('[c]1[c]([OH])[c]([c])[o][c]2[c]([c]([c]1)C(=O)[c]3[c]2[c][c][c][c]3)')

    if not mol.HasSubstructMatch(flavonol_pattern):
        return False, "Not a flavonol structure"


    # If the flavonol pattern matches, it's a flavonol
    return True, "Flavonol structure found"