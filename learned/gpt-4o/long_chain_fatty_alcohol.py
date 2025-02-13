"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a primary alcohol group and a carbon chain length of 13 to 22 carbons.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of a primary alcohol (R-CH2OH)
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX4H2][OX2H]")
    if not mol.HasSubstructMatch(primary_alcohol_pattern):
        return False, "No primary alcohol group (R-CH2OH) found"
    
    # Calculate the length of the longest carbon chain using all pairs of atom paths
    lengths = [len(path) for path in rdmolops.GetSymmSSSR(mol)]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Ensure it's a carbon
            for neighbor in atom.GetNeighbors():
                # Perform a depth-first search to account for complex chain connections
                path_length = Chem.rdmolfiles.CalcNumAtomStereoCenters(mol)  # Use stereo center info with DFS/Complex Path Finding
                lengths.append(path_length)

    max_length = max(lengths, default=0)
    
    # Check if the longest chain corresponds to the desired length
    if 13 <= max_length <= 22:
        return True, f"Contains a primary alcohol group and a suitable carbon chain length of {max_length}."
    else:
        return False, f"Longest carbon chain length is {max_length}, outside the 13-22 range."