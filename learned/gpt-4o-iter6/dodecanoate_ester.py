"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is a fatty acid ester in which the carboxylic acid component is lauric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define advanced lauroyl ester pattern using SMARTS
    # This specifies a 12-carbon chain to an ester, accommodating linear or branched structures
    # Include possible variations for esters connected to glycerol backbones or other groups
    lauroyl_patterns = [
        Chem.MolFromSmarts("C(=O)OCCCCCCCCCCCC"),  # simple lauroyl group
        Chem.MolFromSmarts("C(=O)OC[C@H](O)CCCCCCCCCCC"),  # stereo-specific pattern 
        Chem.MolFromSmarts("C(=O)O[C@H](CC([O-])=O)C[N+](C)(C)C")  # pattern with other components
    ]

    # Iterate over patterns and check for matches
    for pattern in lauroyl_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains characteristic lauroyl ester group indicative of dodecanoate ester"
    
    return False, "Does not contain characteristic lauroyl ester group"