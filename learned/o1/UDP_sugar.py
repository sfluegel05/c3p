"""
Classifies: CHEBI:17297 UDP-sugar
"""
"""
Classifies: UDP-sugar
"""
from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    A UDP-sugar is a pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to a sugar via an anomeric diphosphate linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for UDP moiety (uridine attached to diphosphate)
    udp_smarts = """
    n1(c(=O)[nH]c(=O)c1=O)
    [C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O
    P(=O)(O)OP(=O)(O)O
    """
    udp_smarts = udp_smarts.replace('\n', '')
    udp_pattern = Chem.MolFromSmarts(udp_smarts)
    if udp_pattern is None:
        return False, "Unable to create UDP pattern"

    # Check for UDP moiety
    if not mol.HasSubstructMatch(udp_pattern):
        return False, "UDP moiety not found"

    # Define SMARTS pattern for sugar moiety attached via diphosphate linkage
    # Looking for sugar ring attached to diphosphate oxygen
    sugar_attachment_smarts = "O[P](=O)(O)OP(=O)(O)[O][C]"
    sugar_attachment_pattern = Chem.MolFromSmarts(sugar_attachment_smarts)
    if sugar_attachment_pattern is None:
        return False, "Unable to create sugar attachment pattern"

    # Check for sugar attached via diphosphate linkage
    if not mol.HasSubstructMatch(sugar_attachment_pattern):
        return False, "Sugar moiety not found attached via diphosphate linkage"

    return True, "UDP-sugar identified with UDP moiety linked to sugar via diphosphate"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:46211',
        'name': 'UDP-sugar',
        'definition': 'A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to an unspecified sugar via an anomeric diphosphate linkage.',
        'parents': []
    },
    'config': {}
}