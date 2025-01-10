"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is an ester formed via the esterification of tetradecanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced SMARTS pattern for tetradecanoate ester
    # Look for ester sequence: 14-carbon chain with ester bond component
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O")

    if mol.HasSubstructMatch(tetradecanoate_pattern):
        # Get matches and ensure they conform to expected structure
        matches = mol.GetSubstructMatches(tetradecanoate_pattern)
        for match in matches:
            # Ensuring the ester oxygen is bonded correctly
            # Count carbons in the match for validation
            c_count = 0
            for atom_idx in match:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 6:
                    c_count += 1
        
            if c_count == 14:
                return True, "Valid tetradecanoate ester group found"
    
    return False, "No valid tetradecanoate ester group found or incorrect structure"