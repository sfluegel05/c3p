"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: CHEBI:46617 ribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside based on its SMILES string.
    A ribonucleoside consists of D-ribose connected to a nucleobase via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ribonucleoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Define D-ribose pattern (furanose form with correct stereochemistry)
    ribose_pattern = Chem.MolFromSmarts(
        "[C@H]1(O[C@H](CO)[C@H](O)[C@@H]1O)"
    )
    if not mol.HasSubstructMatch(ribose_pattern):
        return False, "No D-ribose sugar detected"

    # Find the ribose oxygen (O in the ring)
    ribose_matches = mol.GetSubstructMatches(ribose_pattern)
    if not ribose_matches:
        return False, "Ribose structure mismatch"

    # Check glycosidic bond (ribose C1 connected to nucleobase N)
    glycosidic_bond = False
    for match in ribose_matches:
        anomeric_carbon = match[0]  # First atom in SMARTS (C1)
        for neighbor in mol.GetAtomWithIdx(anomeric_carbon).GetNeighbors():
            if neighbor.GetAtomicNum() == 7:  # Connected to nitrogen (nucleobase)
                glycosidic_bond = True
                break
        if glycosidic_bond:
            break

    if not glycosidic_bond:
        return False, "No glycosidic bond to nucleobase"

    # Check for at least one ring aromatic atom (crude nucleobase check)
    has_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())
    if not has_aromatic:
        return False, "No aromatic system (missing nucleobase?)"

    return True, "D-ribose sugar connected to nucleobase via glycosidic bond"