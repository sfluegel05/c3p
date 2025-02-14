"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: CHEBI:31647 7-hydroxyisoflavone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, TautomerEnumerator

def is_7_hydroxyisoflavone(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone is an isoflavone compound with a hydroxy group at the 7-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES and enumerate tautomers
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    tautomers = TautomerEnumerator(mol)
    mol = tautomers.EnumerateIsomericSmiles(False, True)

    # Look for isoflavone backbone pattern
    isoflavone_pattern = Chem.MolFromSmarts("O=C1C=C(O)c2ccccc2OC1")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone backbone found"

    # Look for hydroxy group at 7-position
    hydroxy_7_pattern = Chem.MolFromSmarts("O=C1C=C(O)c2ccc(O)cc2OC1")
    if not mol.HasSubstructMatch(hydroxy_7_pattern):
        return False, "No hydroxy group at the 7-position of the isoflavone backbone"

    # Check for the absence of additional hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("OC")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) > 1:
        return False, "Additional hydroxy groups present"

    # Check if the molecule is a known 7-hydroxyisoflavone
    known_7_hydroxyisoflavones = [
        "CC(C)(O)CCc1cc(ccc1O)-c1coc2cc(O)cc(O)c2c1=O",  # isowigtheone hydrate
        "CC(C)=CCc1c(O)ccc(c1O)-c1coc2cc(O)cc(O)c2c1=O",  # licoisoflavone A
        # Add more known SMILES strings as needed
    ]
    if Chem.MolToSmiles(mol) in known_7_hydroxyisoflavones:
        return True, "Matches a known 7-hydroxyisoflavone structure"

    return True, "Contains an isoflavone backbone with a hydroxy group at the 7-position"