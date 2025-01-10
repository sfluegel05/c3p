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

    # Updated SMARTS for dihydroflavonol core: a generic flavanone core with potential hydroxyl at position 3
    # Accommodate possible stereochemistry and hydroxy position at C3
    core_pattern = Chem.MolFromSmarts("[OH][C@@H]1Cc2c(O)ccc(O)c2/C(=O)[C@H](O1)c3ccccc3")
    
    # Check if the molecule contains the dihydroflavonol core pattern
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Does not match basic dihydroflavonol core structure"

    # Pattern for the potential additional hydroxyl group on the chromanone moiety
    hydroxy_pattern = Chem.MolFromSmarts("O[C@H]1Cc2c(O)ccc(O)c2C1")

    # Check if the molecule has the additional hydroxyl group at the traditional C3 position
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Hydroxyl group missing from expected position"

    # If passes all checks, it's identified as a dihydroflavonol
    return True, "Contains dihydroflavonol core structure with hydroxyl group at position 3"

# Example usage for testing: is_dihydroflavonols("OC1C(Oc2cc(O)cc(O)c2C1=O)c1cc(O)c(O)c(O)c1")