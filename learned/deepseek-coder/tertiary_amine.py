"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:32876 tertiary amine
A tertiary amine is a compound formally derived from ammonia by replacing three hydrogen atoms by hydrocarbyl groups.
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine has a nitrogen atom bonded to three carbon atoms (hydrocarbyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for tertiary amines
    # This matches nitrogen with exactly 3 single bonds to carbon atoms
    # and accounts for both neutral and charged nitrogen
    tertiary_amine_pattern = Chem.MolFromSmarts("[N;X3;+0:1]([C;!$(C=O)])[C;!$(C=O)][C;!$(C=O)] | [N+;X4:1]([C;!$(C=O)])([C;!$(C=O)])[C;!$(C=O)]")
    
    # Check if the molecule matches the pattern
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Contains a nitrogen atom bonded to three carbon atoms (tertiary amine)"
    
    # Check for aromatic nitrogen with three carbon bonds
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            # Get all carbon neighbors
            carbon_neighbors = [neighbor for neighbor in atom.GetNeighbors() 
                              if neighbor.GetAtomicNum() == 6 and 
                              not neighbor.GetIsAromatic() and 
                              not any(bond.GetBondType() == Chem.BondType.DOUBLE 
                                    for bond in mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()))]
            
            # Check if nitrogen has exactly 3 carbon neighbors
            if len(carbon_neighbors) == 3:
                return True, "Contains a nitrogen atom bonded to three carbon atoms (tertiary amine)"

    return False, "No nitrogen atom bonded to three carbon atoms found"