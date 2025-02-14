"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
"""
Classifies: CHEBI:50510 2-oxo monocarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid is a monocarboxylic acid with a 2-oxo substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, "Molecule does not contain exactly one carboxylic acid group"
    
    # Check for presence of 2-oxo group
    oxo_pattern = Chem.MolFromSmarts("C(=O)C")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if len(oxo_matches) == 0:
        return False, "Molecule does not contain a 2-oxo group"
    
    # Check if 2-oxo group and carboxylic acid group are connected
    for oxo_match in oxo_matches:
        oxo_atom_idx = oxo_match[0]
        neighbor_atoms = [mol.GetAtomWithIdx(n).GetSymbol() for n in mol.GetAtomWithIdx(oxo_atom_idx).GetNeighbors()]
        if "O" in neighbor_atoms:
            return True, "Molecule contains a 2-oxo group and a carboxylic acid group connected by a carbon atom"
    
    return False, "2-oxo group and carboxylic acid group are not connected"