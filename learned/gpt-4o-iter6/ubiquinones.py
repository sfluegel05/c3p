"""
Classifies: CHEBI:16389 ubiquinones
"""
from rdkit import Chem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    Ubiquinones are benzoquinones with methoxy groups and a polyisoprenoid side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for general benzoquinone core with flexibility
    benzoquinone_core_pattern = Chem.MolFromSmarts("C1=CC(=O)C=CC(=O)C1")
    if not mol.HasSubstructMatch(benzoquinone_core_pattern):
        return False, "No benzoquinone core moiety found"

    # Check for methoxy groups generally located on the benzoquinone core
    methoxy_match = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() == 1:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and any([n.GetAtomicNum() == 6 for n in neighbor.GetNeighbors()]) == 1:
                    methoxy_match = True
                    break
        if methoxy_match: 
            break
    if not methoxy_match:
        return False, "Required methoxy groups not found on the benzoquinone moiety"

    # Check for extended polyisoprenoid chain. 
    # We will look for repetitive isoprene units.
    polyisoprenoid_pattern = Chem.MolFromSmarts("C=C(C)CCC=C")
    matches = mol.GetSubstructMatches(polyisoprenoid_pattern)
    if len(matches) < 2:
        return False, "No extended polyisoprenoid chain found"
    
    return True, "Matches the characteristics of a ubiquinone: contains a typical benzoquinone moiety with methoxy groups and an extended polyisoprenoid chain"