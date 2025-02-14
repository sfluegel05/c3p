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
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for 2-oxo monocarboxylic acid anion
    # The pattern looks for a carbon with an oxo group in the 2-position, followed by a carboxylate
    oxo_monocarboxylic_pattern = Chem.MolFromSmarts("C(=O)[C;R0][CH2,C][C;R0](=O)[O-]")
    
    # Check if the pattern matches within the molecule
    if mol.HasSubstructMatch(oxo_monocarboxylic_pattern):
        return True, "Contains a 2-oxo group adjacent to a carboxylate group in the correct position"
    
    # Additional verification for complex structures
    # Check carbon positions and connectivity
    for atom in mol.GetAtoms():
        # Check for oxo group at 2-position only if it's part of a C-C-C=O sequence
        if atom.GetAtomicNum() == 6 and len(atom.GetNeighbors()) >= 3:
            oxy_groups = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()).GetBondType().name == 'DOUBLE']
            carboxylate_groups = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and any([n.GetSymbol() == 'O' for n in nbr.GetNeighbors()])]

            if len(oxy_groups) == 1 and len(carboxylate_groups) == 1:
                return True, "Contains a 2-oxo group adjacent to a carboxylate group in a complex position"
    
    return False, "Does not contain a 2-oxo group adjacent to a carboxylate group"