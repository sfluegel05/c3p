"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    
    A 2-hydroxydicarboxylic acid contains two carboxylic acid groups and a hydroxy 
    group on the carbon atom at position alpha to the carboxy group.

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

    # Define the carboxylic acid pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    if len(carboxylic_acid_matches) < 2:
        return False, "Found fewer than two separate carboxylic acid groups"
    
    # Analyze potential alpha-carbon with hydroxyl attachment for each carboxylic group
    for idx1, carbox1 in enumerate(carboxylic_acid_matches):
        for idx2, carbox2 in enumerate(carboxylic_acid_matches):
            if idx1 >= idx2:  # Avoid duplicate pair analysis and self-pair analysis
                continue
            
            carbox1_carbon = carbox1[0]
            carbox2_carbon = carbox2[0]

            # Check adjacent carbons (potential alpha carbons)
            candidate_alpha_carbons = set(atom.GetIdx() for atom in mol.GetAtomWithIdx(carbox1_carbon).GetNeighbors()) & \
                                      set(atom.GetIdx() for atom in mol.GetAtomWithIdx(carbox2_carbon).GetNeighbors())
            
            for alpha_carbon in candidate_alpha_carbons:
                # Check for hydroxyl group on the candidate alpha carbon
                for neighbor in mol.GetAtomWithIdx(alpha_carbon).GetNeighbors():
                    if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(alpha_carbon, neighbor.GetIdx()).GetBondTypeAsDouble() == 1.0:
                        return True, "Contains two carboxylic acid groups and a hydroxy group on the alpha carbon"

    return False, "No suitable 2-hydroxydicarboxylic acid pattern found"

# Example test cases
print(is_2_hydroxydicarboxylic_acid("CC(C(O)=O)C(C)(O)C(O)=O"))  # 2,3-dimethylmalic acid