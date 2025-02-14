"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid (SCFA) based on its SMILES string.
    An SCFA is defined as an aliphatic monocarboxylic acid with a chain length of less than C6.
    Acceptable substituents are limited to hydroxy and keto groups (oxygen-containing).
    
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
    
    # Check for exactly one carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 1:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 1"
    
    # Ensure molecule is aliphatic (no aromatic atoms)
    if any(atom.GetIsAromatic() for atom in mol.GetAtoms()):
        return False, "Molecule contains aromatic atoms"
    
    # Ensure all atoms are C, H, or O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:
            return False, f"Contains atom other than C, H, O: {atom.GetSymbol()}"
    
    # Function to find the longest carbon chain recursively
    def longest_chain_length(mol, atom_idx, visited):
        """
        Recursively find the longest carbon chain starting from atom_idx.
        """
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() != 6:  # Not a carbon atom
            return 0
        max_length = 0
        visited.add(atom_idx)
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in visited and neighbor.GetAtomicNum() == 6:
                length = longest_chain_length(mol, n_idx, visited.copy())
                if length > max_length:
                    max_length = length
        return 1 + max_length  # Include current atom
    
    # Get indices of carbon atoms
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return False, "No carbon atoms found in molecule"
    
    # Find the longest carbon chain
    max_chain_length = 0
    for idx in carbon_atoms:
        chain_length = longest_chain_length(mol, idx, set())
        if chain_length > max_chain_length:
            max_chain_length = chain_length
    
    # Check if chain length is less than 6
    if max_chain_length >= 6:
        return False, f"Longest carbon chain length is {max_chain_length}, which is not less than 6"
    
    return True, "Molecule is an aliphatic monocarboxylic acid with chain length less than C6"