"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are flavan-3-ol compounds characterized by a flavan skeleton and hydroxylation patterns.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """

    # Convert SMILES string to rdkit Molecule object
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for catechins, focusing on flavan-3-ol core
    # The pattern matches the structural motif: 2-phenyl-3,4-dihydro-2H-chromene (flavan-3-ol core)
    # The pattern is broader to accommodate various substituted catechins and stereoisomers
    catechin_pattern = Chem.MolFromSmarts("C1COc2cccc(O)c2[C@H]1Oc3ccc(O)c(O)c3")

    if not mol.HasSubstructMatch(catechin_pattern):
        return False, "No flavan-3-ol skeleton match found"

    # Further refinement of pattern matching can include specific hydroxylation if necessary
    return True, "Matches flavan-3-ol skeleton of catechins"

# The revised code now employs a broader SMARTS pattern to account for catechin structures, accommodating potential isomers and modifications typical of catechin derivatives.