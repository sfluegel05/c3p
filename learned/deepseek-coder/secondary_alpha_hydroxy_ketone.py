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

    # Define a more flexible SMARTS pattern for secondary alpha-hydroxy ketone
    # This pattern looks for:
    # - A carbon with a hydroxyl group (C-OH)
    # - Adjacent to a carbonyl carbon (C=O)
    # - The hydroxyl-bearing carbon is bonded to at least one non-hydrogen atom (organyl group)
    pattern = Chem.MolFromSmarts("[C;H0,H1]([OH])[C](=O)[C;!H0]")

    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(pattern):
        # Verify that the hydroxyl-bearing carbon is a secondary carbon
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            hydroxyl_carbon_idx = match[0]
            hydroxyl_carbon = mol.GetAtomWithIdx(hydroxyl_carbon_idx)
            
            # Count non-hydrogen neighbors (should be at least 2: OH and organyl group)
            non_h_neighbors = [n for n in hydroxyl_carbon.GetNeighbors() if n.GetAtomicNum() != 1]
            if len(non_h_neighbors) >= 2:
                # Check that one neighbor is OH and another is an organyl group
                has_oh = any(n.GetAtomicNum() == 8 and n.GetTotalNumHs() == 1 for n in non_h_neighbors)
                has_organyl = any(n.GetAtomicNum() != 8 for n in non_h_neighbors)
                
                if has_oh and has_organyl:
                    return True, "Contains a secondary alpha-hydroxy ketone (acyloin) structure"
        
        return False, "Hydroxyl-bearing carbon does not meet secondary alpha-hydroxy ketone criteria"
    else:
        return False, "Does not contain a secondary alpha-hydroxy ketone structure"