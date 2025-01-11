"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:<id> dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is an ester where the carboxylic acid component is lauric acid (dodecanoic acid).

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

    # Define lauric acid molecule for comparison
    lauric_acid = Chem.MolFromSmiles('CCCCCCCCCCC(=O)O')  # Lauric acid SMILES

    # Define ester hydrolysis reaction to get carboxylic acids
    rxn = AllChem.ReactionFromSmarts('[C:1](=O)[O][C:2]>>[C:1](=O)[O]')  # Simplified ester hydrolysis

    # Run the reaction on the molecule
    products = rxn.RunReactants((mol,))

    # Keep track of whether lauric acid is found
    found_lauric_acid = False

    # Iterate over reaction products
    for product_set in products:
        acid = product_set[0]  # Get the acid fragment
        # Sanitize the fragment
        Chem.SanitizeMol(acid, catchErrors=True)
        # Check if the fragment matches lauric acid
        if acid.HasSubstructMatch(lauric_acid):
            found_lauric_acid = True
            break  # No need to check further

    if found_lauric_acid:
        return True, "Contains lauric acid ester group"
    else:
        return False, "No lauric acid ester groups found"