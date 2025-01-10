"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins are flavan-3-ol compounds with a flavan skeleton and hydroxylation patterns.
    
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

    # Define a more specific SMARTS pattern for catechins considering hydroxyl and stereochemistry
    # Pattern: aromatic A and B rings, central C ring with hydroxyls, stereochemistry at position 3
    flavan_3_ol_pattern = Chem.MolFromSmarts("c1c(O)cc(O)cc1C2(CO[C@@H]2c3ccc(O)c(O)c3)")

    if not mol.HasSubstructMatch(flavan_3_ol_pattern):
        return False, "No flavan-3-ol skeleton match found"

    # Further checks for the presence of typical catechin hydroxylation patterns can be added here
    
    return True, "Matches flavan-3-ol skeleton of catechins"

# Note: This pattern checks for the flavan-3-ol structure, considering typical stereochemistry and hydroxylation.
# Testing on a variety of examples can help refine it further to cover all catechin derivatives.