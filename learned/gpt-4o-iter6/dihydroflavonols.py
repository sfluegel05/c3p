"""
Classifies: CHEBI:48039 dihydroflavonols
"""
from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    Dihydroflavonols are characterized by a hydroxyflavanone structure with a hydroxy group at position 3.

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

    # Define a more flexible SMARTS pattern for the dihydroflavonol core
    # Allow for possible variance in stereochemistry and flexible position of hydroxyl groups
    core_pattern = Chem.MolFromSmarts("O[C@@H]1Cc2cc(O)ccc2/C(=O)[C@H](O1)c1cc(O)ccc1")
    
    # Check if the molecule contains the dihydroflavonol core pattern
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Does not match generalized dihydroflavonol core structure"

    # Pattern for the presence of hydroxy groups; focus less on specific atom indexing to match class characteristics
    hydroxy_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("O"))
    if len(hydroxy_groups) < 3:
        return False, "Less than three hydroxy groups detected, typical for dihydroflavonols"

    # If passes all checks, it's identified as a dihydroflavonol
    return True, "Contains dihydroflavonol core structure with sufficient hydroxylation"

# Example for testing the function
# is_dihydroflavonols("OC1C(Oc2cc(O)cc(O)c2C1=O)c1cc(O)c(O)c(O)c1")