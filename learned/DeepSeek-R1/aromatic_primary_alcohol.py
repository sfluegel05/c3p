"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: CHEBI:134307 aromatic primary alcohol
"""
from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has a hydroxyl group (-OH) attached to a primary carbon (CH2)
    that is directly bonded to an aromatic carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate through all oxygen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            # Check if it's a hydroxyl group (-OH)
            if atom.GetDegree() == 1 and atom.GetTotalNumHs() == 1:
                # Get the connected carbon (primary alcohol carbon)
                carbon = atom.GetNeighbors()[0]
                # Check if it's a primary alcohol (CH2OH)
                if carbon.GetDegree() == 2 and carbon.GetAtomicNum() == 6:
                    # Get the adjacent carbon (attached to the CH2 group)
                    adjacent_carbons = [n for n in carbon.GetNeighbors() if n != atom]
                    if len(adjacent_carbons) != 1:
                        continue  # Shouldn't happen due to degree check
                    adjacent_carbon = adjacent_carbons[0]
                    
                    # Check if the adjacent carbon is aromatic
                    if adjacent_carbon.GetIsAromatic():
                        return True, "Primary alcohol group attached to aromatic carbon"
    
    return False, "No primary alcohol group attached to aromatic carbon found"