"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:36299 dodecanoate ester
Any fatty acid ester in which the carboxylic acid component is lauric acid (dodecanoic acid).
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.

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

    # Look for dodecanoate fragment (CCCCCCCCCCCCC(=O)O)
    dodecanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCCCC(=O)O")
    dodecanoate_matches = mol.GetSubstructMatches(dodecanoate_pattern)
    if not dodecanoate_matches:
        return False, "No dodecanoate fragment found"

    # Look for ester bond (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester bond found"

    # Check that dodecanoate is the carboxylic acid component
    for match in dodecanoate_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == "O":
                neighbors = atom.GetNeighbors()
                if len(neighbors) == 1:
                    neighbor = neighbors[0]
                    if neighbor.GetSymbol() == "C":
                        neighbor_neighbors = neighbor.GetNeighbors()
                        if len(neighbor_neighbors) == 3:
                            for nn in neighbor_neighbors:
                                if nn.GetSymbol() == "O" and nn.GetFormalCharge() == 0:
                                    return True, "Contains dodecanoate as the carboxylic acid component of an ester"

    return False, "Dodecanoate not found as the carboxylic acid component of an ester"