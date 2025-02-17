"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: Epoxide
Definition: Any cyclic ether in which the oxygen atom forms part of a 3-membered ring.
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if the molecule represented by the SMILES string is an epoxide.
    An epoxide is defined as any cyclic ether in which the oxygen atom forms part of a 3-membered ring
    (i.e. a ring containing exactly one oxygen and two carbons).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as an epoxide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the provided SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a 3-membered cyclic ether:
    # The pattern looks for a three-membered ring ([R3]) that contains two carbon atoms and one oxygen atom.
    epoxide_pattern = Chem.MolFromSmarts("[C;R3]-[O;R3]-[C;R3]")
    if epoxide_pattern is None:
        return False, "Failed to generate SMARTS pattern for epoxide"

    # Find substructure matches: this will return a tuple for every match found
    matches = mol.GetSubstructMatches(epoxide_pattern)
    if matches:
        return True, f"Found epoxide ring(s): {len(matches)} occurrence(s) of a 3-membered cyclic ether"
    else:
        return False, "No epoxide ring (3-membered cyclic ether with one oxygen and two carbons) found"