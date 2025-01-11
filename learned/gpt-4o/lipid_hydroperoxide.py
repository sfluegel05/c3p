"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is a lipid molecule with one or more hydroperoxy (-OOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroperoxy group pattern (-O-O)
    hydroperoxy_pattern = Chem.MolFromSmarts("OO")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy group (-OOH) found"

    # Check the attachment of the hydroperoxy group to long hydrocarbon chain
    # Typically expect a hydrocarbon chain like -[CH2]-[CH]=[CH]- with -OOH attached
    carbon_chain_pattern = Chem.MolFromSmarts("[CH2,CH]=[CH,CH2]-[CH,CH2]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Hydroperoxy group not attached to a suitable carbon chain indicative of lipids"

    return True, "Contains hydroperoxy groups attached to a lipid backbone"