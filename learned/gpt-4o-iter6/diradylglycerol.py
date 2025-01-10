"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol has a glycerol backbone with exactly two substituent groups (acyl, alkyl, or alk-1-enyl).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved glycerol backbone check, allowing for stereo
    glycerol_pattern = Chem.MolFromSmarts("C[C@H](O)CO")  # Recognize one of the chiral carbons in glycerol
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone pattern not matched"

    # Recognize ester link (acyl group: R-C=O-), ether link (alkyl: R-O-), or alk-1-enyl pattern
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[O][CH2]")
    ether_pattern = Chem.MolFromSmarts("[OX2][CX4]")
    alk1enyl_pattern = Chem.MolFromSmarts("[CX3](=[CH2])[CH](O)[CH2]")

    # Finding substituents by these patterns
    # Limit to substituents that connect via ester or ether bonds
    acyl_matches = mol.GetSubstructMatches(ester_pattern)
    alkyl_matches = mol.GetSubstructMatches(ether_pattern)
    alk1enyl_matches = mol.GetSubstructMatches(alk1enyl_pattern)

    total_unique_substituents = (
        len(set([match[0] for match in acyl_matches])) +
        len(set([match[0] for match in alkyl_matches])) +
        len(set([match[0] for match in alk1enyl_matches]))
    )

    # Ensure there are exactly two substituents
    if total_unique_substituents != 2:
        return False, f"Expected exactly 2 substituent groups, found {total_unique_substituents}"

    return True, "Valid diradylglycerol with two substituents connected via ester, ether, or alk-1-enyl bonds"