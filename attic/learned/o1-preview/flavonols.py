"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: flavonols
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is defined as any hydroxyflavone where the hydrogen at position 3
    of the heterocyclic ring is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the flavonol core (3-hydroxyflavone)
    # SMILES for 3-hydroxyflavone
    flavonol_smiles = 'O=C1c2ccccc2Oc2cc(O)ccc12'
    flavonol_core = Chem.MolFromSmiles(flavonol_smiles)
    if flavonol_core is None:
        return False, "Unable to construct flavonol core structure"

    # Check for flavonol core substructure
    if not mol.HasSubstructMatch(flavonol_core):
        return False, "Flavonol core structure not found"

    return True, "Contains flavonol core structure with 3-hydroxy group"

__metadata__ = {
    'chemical_class': {
        'name': 'flavonols',
        'definition': 'Any hydroxyflavone in which the ring hydrogen at position 3 of the heterocyclic ring is replaced by a hydroxy group.',
    },
    'config': {
        # Configuration and performance metrics can be added here if needed
    },
    # Additional metadata can be added here
}