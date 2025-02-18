"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: CHEBI:36699 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is formed by condensation of a carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for exactly one ester group
    ester_pattern = Chem.MolFromSmarts("[O;X2][C;X3](=[O;X1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Expected 1 ester group, found {len(ester_matches)}"

    # Get the oxygen and adjacent carbon in sterol part
    oxygen_idx, carbonyl_idx = ester_matches[0][0], ester_matches[0][1]
    oxygen = mol.GetAtomWithIdx(oxygen_idx)
    sterol_carbon = next(n for n in oxygen.GetNeighbors() if n.GetIdx() != carbonyl_idx)

    # Check if sterol carbon is part of a steroid nucleus (tetracyclic system)
    steroid_smarts = Chem.MolFromSmarts(
        "[C]1CC[C@]2[C@@H]3CC=C4C[C@H](CC[C@]4(C3)CC[C@]12C)"
    )
    if not mol.HasSubstructMatch(steroid_smarts):
        return False, "No steroid nucleus detected"

    # Verify sterol carbon is in the matched steroid structure
    steroid_matches = mol.GetSubstructMatches(steroid_smarts)
    steroid_atoms = {idx for match in steroid_matches for idx in match}
    if sterol_carbon.GetIdx() not in steroid_atoms:
        return False, "Ester not attached to steroid nucleus"

    # Check acid chain: carbonyl must connect to a carbon chain
    carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_idx)
    acid_chain = any(n.GetAtomicNum() == 6 for n in carbonyl_carbon.GetNeighbors() if n.GetIdx() != oxygen_idx)
    if not acid_chain:
        return False, "No acid chain attached to ester"

    return True, "Steroid nucleus with ester-linked acid chain"