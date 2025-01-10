"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is an ester formed from tetradecanoic acid (14-carbon chain).

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

    # Define an enhanced SMARTS pattern for a tetradecanoate ester
    # Looking for ester group with a 14-carbon chain specifically
    tetradecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC(CC)=O")  # Sequence for 14-carbons ending in carboxylic acid esters
   
    if mol.HasSubstructMatch(tetradecanoate_pattern):
        # Get matches and ensure they conform to expected carbon count
        matches = mol.GetSubstructMatches(tetradecanoate_pattern)
        for match in matches:
            # The ends of the matched SMARTS should indicate attachment to ester oxygen
            c_count = 0
            for atom_idx in match:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 6:
                    c_count += 1

            if c_count == 14:  # Ensuring the exact carbon count matching myristic acid
                return True, "Valid tetradecanoate ester group found"
    
    return False, "No valid tetradecanoate ester group found or incorrect carbon chain"