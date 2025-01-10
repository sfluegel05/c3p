"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is a monoglyceride in which the acyl substituent is located at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the monoglyceride pattern with acyl group at position 1
    monoglyceride_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0][CH2][CH](O)[CH2]O")
    if not mol.HasSubstructMatch(monoglyceride_pattern):
        return False, "Molecule does not match 1-monoglyceride pattern"

    # Ensure there is exactly one ester group
    ester_pattern = Chem.MolFromSmarts("[$([CX3](=O)[OX2H1]),$([CX3](=O)[O–][#6])]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check that positions 2 and 3 have free hydroxyl groups
    glycerol_oh_pattern = Chem.MolFromSmarts("[CH2][CH](O)[CH2]O")
    if not mol.HasSubstructMatch(glycerol_oh_pattern):
        return False, "Positions 2 and/or 3 do not have free hydroxyl groups"

    # Check that there are no additional acyl chains or ester groups
    # Count total number of acyl chains (ester and amide bonds)
    acyl_chain_count = len(ester_matches)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    acyl_chain_count += len(amide_matches)
    if acyl_chain_count != 1:
        return False, f"Molecule has {acyl_chain_count} acyl chains, expected 1"

    # Check molecular weight - monoglycerides typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for a 1-monoglyceride"

    return True, "Molecule is a 1-monoglyceride with acyl group at position 1 and free hydroxyls at positions 2 and 3"