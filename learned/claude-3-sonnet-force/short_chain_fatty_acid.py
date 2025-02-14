"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: CHEBI:35834 short-chain fatty acid
A short-chain fatty acid is an aliphatic monocarboxylic acid with a chain length of less than C6.
If any non-hydrocarbon substituent is present, the compound is not normally regarded as a short-chain fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Convert to SMARTS pattern
    smarts = Chem.MolToSmarts(mol)
    
    # Check if SMARTS matches pattern for short-chain fatty acid
    pattern = r"C(=O)O[C@,C@@,C]([C,c])[C,c]([C,c])([C,c])=O"
    if not Chem.MolFromSmarts(pattern).HasSubstructMatch(mol):
        return False, "Molecule does not match pattern for short-chain fatty acid"
    
    # Check for non-hydrocarbon substituents (allow hydroxy groups)
    non_h_substituents = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != 'H' and atom.GetSymbol() != 'C' and atom.GetSymbol() != 'O':
            non_h_substituents = True
            break
        elif atom.GetSymbol() == 'C' or atom.GetSymbol() == 'O':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() not in ['H', 'C', 'O']:
                    non_h_substituents = True
                    break
    if non_h_substituents:
        return False, "Molecule contains non-hydrocarbon substituents"
    
    return True, "Aliphatic monocarboxylic acid with chain length < C6 and no non-hydrocarbon substituents (except hydroxy groups)"