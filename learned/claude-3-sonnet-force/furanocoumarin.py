"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: CHEBI:27727 furanocoumarin

A furanocoumarin is defined as any furochromene that consists of a furan ring fused with a coumarin.
The fusion may occur in different ways to give several isomers.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_furanocoumarin(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for the fused furanocoumarin system
    furanocoumarin_pattern = Chem.MolFromSmarts("c1cc2c3c(c1)oc1ccccc1c3oc2")

    if mol.HasSubstructMatch(furanocoumarin_pattern):
        # Additional checks for furanocoumarin characteristics (optional)
        # ...

        return True, "Contains a fused furanocoumarin system"

    return False, "Does not match the structural characteristics of a furanocoumarin"