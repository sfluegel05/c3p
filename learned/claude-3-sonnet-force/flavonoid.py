"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: CHEBI:35456 flavonoid

A flavonoid is any member of the 'superclass' flavonoids whose skeleton is based on 1-benzopyran
with an aryl substituent at position 2. The term was originally restricted to natural products,
but is now also used to describe semi-synthetic and fully synthetic compounds.
"""

from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """

    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Define flavonoid scaffold pattern
        scaffold_pattern = Chem.MolFromSmarts("[c&r1]1[c&r2]c2c(cc1)oc(c2)-[c&r3]")

        # Check if scaffold is present
        if mol.HasSubstructMatch(scaffold_pattern):
            return True, "Contains flavonoid scaffold"
        else:
            return False, "Flavonoid scaffold not present"

    except Exception as e:
        return False, f"Error: {str(e)}"