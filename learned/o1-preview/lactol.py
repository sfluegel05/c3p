"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: lactol
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group to an aldehyde or ketone.
    They are 1-oxacycloalkan-2-ols or unsaturated analogues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for lactol
    # Pattern explanation:
    # - Ring structure containing O (ether oxygen)
    # - Hydroxyl group [OH] attached to adjacent carbon
    # - Adjacent carbon to oxygen in the ring
    # - Exclude ethers (no carbonyl involved)
    lactol_pattern = Chem.MolFromSmarts("""
    [#6;R][O;R][#6;R]([OH])[#6;R]
    """)
    if lactol_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for lactol substructure
    matches = mol.GetSubstructMatches(lactol_pattern)
    if matches:
        return True, "Contains lactol moiety (cyclic hemiacetal)"
    else:
        return False, "No lactol moiety found"

__metadata__ = {
    'chemical_class': {
        'name': 'lactol',
        'definition': 'Cyclic hemiacetals formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group. They are thus 1-oxacycloalkan-2-ols or unsaturated analogues.'
    },
}