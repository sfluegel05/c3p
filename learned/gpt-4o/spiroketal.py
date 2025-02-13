"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is defined by a ketal carbon that is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a spiroketal, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a spiroketal: a spirocarbon with two alkoxy groups
    spiroketal_pattern = Chem.MolFromSmarts("[C]1([O])([O])-[R]2-[R]1")

    # Check for spiroketal substructure in the molecule
    if not mol.HasSubstructMatch(spiroketal_pattern):
        return False, "No spiroketal structure found"

    # Verify that two rings share the ketal carbon
    matches = mol.GetSubstructMatches(spiroketal_pattern)
    for match in matches:
        ketal_carbon = match[0]
        neighbor_idxs = [a.GetIdx() for a in mol.GetAtomWithIdx(ketal_carbon).GetNeighbors()]
        rings = [ring for ring in Chem.GetSymmSSSR(mol) if ketal_carbon in ring]
        if len(rings) != 2 or not all(set(neighbor_idxs).intersection(ring) for ring in rings):
            return False, "No valid spiroketal: ketal carbon shared by two rings not found"

    return True, "Contains ketal carbon common to two rings (spiroketal structure)"