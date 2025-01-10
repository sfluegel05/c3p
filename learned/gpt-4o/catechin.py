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

    # Define a broader SMARTS pattern for catechins with a focus on the flavan-3-ol core
    # Allowing flexibility in hydroxylation and aromatic substitution patterns
    catechin_pattern = Chem.MolFromSmarts("c1(c2ccc(OC)c(O)c2)cc3CC(O)C(O)Oc3c1")

    if not mol.HasSubstructMatch(catechin_pattern):
        return False, "No flavan-3-ol skeleton match found"

    # Additional checks for common functional group attachments can be added here if necessary

    return True, "Matches flavan-3-ol skeleton of catechins"

# Note: SMARTS pattern is updated to focus on structure rather than specific stereochemistry,
# to accommodate the stereochemical diversity in catechins effectively.