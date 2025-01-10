"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    
    A 2-hydroxydicarboxylic acid contains two carboxylic acid groups and a hydroxy 
    group on the carbon atom at position alpha to at least one of the carboxy groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxylic acid SMARTS pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)

    # We need to find at least two separate carboxylic acid groups
    if len(carboxylic_acid_matches) < 2:
        return False, "Found fewer than two separate carboxylic acid groups"

    # Hydroxy group pattern 
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4H1,R0]=[OH1]") 
    # Specify precisely for the attachment to non-ring atoms
     
    for carbox in carboxylic_acid_matches:
        carbox_carbon = carbox[0]
        # Find alpha-carbon candidate neighbors from carboxyl carbon that host -OH
        for atom in mol.GetAtomWithIdx(carbox_carbon).GetNeighbors():
            if atom.GetAtomicNum() == 6:  # Ensure it is bonded to a carbon
                alpha_carbon = atom.GetIdx()
                # Check if alpha-carbon is bound to a hydroxy group
                if mol.GetAtomWithIdx(alpha_carbon).GetSmarts() == "[CX4H1]" and any(neigh.GetAtomicNum() == 8 for neigh in mol.GetAtomWithIdx(alpha_carbon).GetNeighbors()):
                    return True, "Contains required carboxylic groups and hydroxy group on the alpha carbon"

    return False, "No suitable 2-hydroxydicarboxylic acid pattern found"