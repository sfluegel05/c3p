"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    Dihydroflavonols have a hydroxyflavanone structure with a hydroxy group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for dihydroflavonol core: 2-phenyl-3-chromanol with hydroxyl at C3 (stereochemistry agnostic)
    core_pattern = Chem.MolFromSmarts("Oc1c(C(=O)[C@H]([C@@H]1O)c2ccccc2)c3ccccc3")

    # Check if the molecule contains the dihydroflavonol core pattern
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Does not match basic dihydroflavonol core structure"

    # Additional hydroxyl group checks are inherently part of the above SMARTS pattern

    # If passes all checks, it's identified as a dihydroflavonol
    return True, "Contains dihydroflavonol core structure with hydroxyl group at position 3"

# Example usage for testing: is_dihydroflavonols("OC1C(Oc2cc(O)cc(O)c2C1=O)c1cc(O)c(O)c(O)c1")