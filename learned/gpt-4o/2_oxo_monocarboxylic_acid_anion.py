"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: 2-Oxo Monocarboxylic Acid Anion
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for 2-oxo monocarboxylic acid anion
    # Looking for a 2-position oxo group and adjacent carboxylate
    oxo_monocarboxylic_pattern = Chem.MolFromSmarts("C([C;R0](=O)[O-])C(=O)")
    
    # Check if the pattern matches within the molecule
    if mol.HasSubstructMatch(oxo_monocarboxylic_pattern):
        return True, "Contains a 2-oxo group adjacent to a carboxylate group in the correct position"
    
    # Additional checks for more complex structure confirmations
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Check for carbon atoms
            # Check for adjacent double bonded oxygen indicating oxo group
            oxo_neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() 
                             if nbr.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()).GetBondType().name == 'DOUBLE']
            
            # Check for carboxylate groups in neighboring carbons
            carboxylate_neighbors = []
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:
                    if any(n.GetAtomicNum() == 8 for n in nbr.GetNeighbors()):
                        carboxylate_neighbors.append(nbr.GetIdx())

            if oxo_neighbors and len(oxo_neighbors) == 1 and carboxylate_neighbors:
                return True, "Contains a 2-oxo group adjacent to a carboxylate group"
    
    return False, "Does not contain a 2-oxo group adjacent to a carboxylate group"