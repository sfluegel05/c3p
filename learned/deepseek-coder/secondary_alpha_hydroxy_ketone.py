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
    # The pattern is more flexible to account for ring systems and complex organyl groups.
    pattern = Chem.MolFromSmarts("[C;H1]([OH])[C](=O)[C]")

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(pattern):
        # Further verify that the hydroxyl-bearing carbon is a secondary carbon
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            hydroxyl_carbon_idx = match[0]
            hydroxyl_carbon = mol.GetAtomWithIdx(hydroxyl_carbon_idx)
            if hydroxyl_carbon.GetDegree() == 3:  # Secondary carbon has 3 bonds (H, OH, and organyl group)
                return True, "Contains a secondary alpha-hydroxy ketone (acyloin) structure"
        
        return False, "Hydroxyl-bearing carbon is not a secondary carbon"
    else:
        return False, "Does not contain a secondary alpha-hydroxy ketone structure"