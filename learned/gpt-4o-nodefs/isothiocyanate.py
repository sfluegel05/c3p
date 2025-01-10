"""
Classifies: CHEBI:52221 isothiocyanate
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    Isothiocyanates ideally have an accessible, reactive functional group N=C=S.

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an isothiocyanate, otherwise False
        str: Reason for classification
    """
    # Parse the SMILES string to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Primary SMARTS pattern for isothiocyanate group: N=C=S
    isothiocyanate_pattern = Chem.MolFromSmarts("N=C=S")
    if mol.HasSubstructMatch(isothiocyanate_pattern):
        # Check whether the N=C=S group is attached to a simple carbon chain
        simple_attachment_pattern = Chem.MolFromSmarts("C(C)N=C=S")
        if mol.HasSubstructMatch(simple_attachment_pattern):
            return True, "Contains accessible isothiocyanate functional group with simple attachments"
        else:
            # Return with an explanation if more complex structures are identified
            return True, "Contains isothiocyanate functional group, though in complex structural context"
    else:
        return False, "Does not contain isothiocyanate functional group (N=C=S)"