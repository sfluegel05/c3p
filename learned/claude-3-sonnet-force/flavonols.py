"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: CHEBI:17773 flavonols
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a hydroxyflavone with a hydroxy group at position 3 of the heterocyclic ring.

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

    # Check if molecule is a flavonoid (has benzopyran backbone)
    flavonoid_pattern = Chem.MolFromSmarts("[c1]2[c]([c]3[c]([c]([c]2[o]1)[O])[O])-[c]([c]([c]3=O)[O])")
    if not mol.HasSubstructMatch(flavonoid_pattern):
        return False, "Not a flavonoid (missing benzopyran backbone)"

    # Check for hydroxy group at position 3
    position_3_oh = Chem.MolFromSmarts("[c1]2[c]([c]([c]([o]2)[OH])[O])-[c]([c]([c]1=O)[O])")
    if not mol.HasSubstructMatch(position_3_oh):
        return False, "No hydroxy group at position 3 of heterocyclic ring"

    return True, "Contains a hydroxyflavone backbone with a hydroxy group at position 3"