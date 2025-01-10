"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: CHEBI:XXXXX secondary alpha-hydroxy ketone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has a hydroxyl group (-OH) and a carbonyl group (C=O) on adjacent carbons,
    where the carbon bearing the hydroxyl group is also bonded to one hydrogen and one organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a secondary alpha-hydroxy ketone
    # The pattern looks for a carbon with a hydroxyl group (C-OH) adjacent to a carbonyl carbon (C=O),
    # where the hydroxyl-bearing carbon is also bonded to one hydrogen and one organyl group.
    pattern = Chem.MolFromSmarts("[CX4H][C](=O)[C]([OH])[CX4]")

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(pattern):
        return True, "Contains a secondary alpha-hydroxy ketone (acyloin) structure"
    else:
        return False, "Does not contain a secondary alpha-hydroxy ketone structure"