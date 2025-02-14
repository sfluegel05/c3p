"""
Classifies: CHEBI:35681 secondary alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule is a secondary alcohol based on its SMILES string.
    A secondary alcohol is defined as a compound in which a hydroxy group (-OH) is attached to a saturated carbon atom
    which has two other carbon atoms attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a secondary alcohol
    secondary_alcohol_pattern = Chem.MolFromSmarts("[C](-[OX2])(-[C])(-[C])")
    
    # Find all occurrences of the pattern in the molecule
    matches = mol.GetSubstructMatches(secondary_alcohol_pattern)

    if not matches:
        return False, "No secondary alcohol substructure found"
    
    for match in matches:
        # Get the central carbon atom index
        central_carbon_idx = match[0]
        central_carbon = mol.GetAtomWithIdx(central_carbon_idx)

        # Get neighbors of the central carbon
        neighbors = [atom for atom in central_carbon.GetNeighbors()]
        carbon_neighbors = 0
        
        # Check neighbors.
        for neighbor in neighbors:
          if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != match[1]: #avoid the oxygen atom
            carbon_neighbors +=1
        if carbon_neighbors == 2:    
          return True, "Molecule contains a secondary alcohol group."
    
    return False, "No secondary alcohol substructure found"