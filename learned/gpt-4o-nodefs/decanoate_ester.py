"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester specifically incorporates a 10-carbon alkyl chain linked via an ester group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a specific decanoate structure using broader pattern recognition 
    # This pattern should account for the ester linkage rO-C(=O)-CCCCCCCCCC which represents the decanoate framework
    # and accommodate permutations and ester group positions in larger structures
    
    # ether and carboxylate variations
    decanoate_patterns = [
        Chem.MolFromSmarts("CCCCCCCCCCCC(=O)O"),  # standard decanoate ester backbone
        Chem.MolFromSmarts("C(=O)OCCCCCCCCCC"),   # reversed ester representation
        Chem.MolFromSmarts("CCCCCCCCCC(=O)O[R]")  # ester linkage as part of larger structures
    ]

    # Attempt to match any of the known patterns
    for pattern in decanoate_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains valid decanoate ester structure"

    return False, "No decanoate ester pattern found"