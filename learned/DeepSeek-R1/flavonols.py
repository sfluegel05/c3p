"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: CHEBI:flavonols
Flavonols are hydroxyflavones with a hydroxyl group at position 3 of the heterocyclic ring.
"""
from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol based on its SMILES string.
    A flavonol is a hydroxyflavone with a hydroxyl group at position 3 of the heterocyclic (C) ring.

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

    # Define SMARTS pattern for flavonol core structure:
    # Benzopyran-4-one (flavone backbone) with hydroxyl at position 3 of the C ring
    flavonol_pattern = Chem.MolFromSmarts("c1ccc2c(c1)oc(=O)c(O)cc2")
    
    # Check for the presence of the core structure
    if mol.HasSubstructMatch(flavonol_pattern):
        return True, "Contains flavone backbone with hydroxyl group at position 3 of the heterocyclic ring"
    else:
        return False, "Does not contain required flavone structure with hydroxyl at position 3"