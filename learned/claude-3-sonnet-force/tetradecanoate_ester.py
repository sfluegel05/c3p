"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    A tetradecanoate ester is a fatty acid ester obtained by condensation of the
    carboxy group of tetradecanoic acid (myristic acid) with a hydroxy group of
    an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetradecanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for tetradecanoic acid fragment
    tetradecanoic_acid = Chem.MolFromSmarts("CCCCCCCCCCCCCC(=O)O")
    tetradecanoic_acid_matches = mol.GetSubstructMatches(tetradecanoic_acid)
    if not tetradecanoic_acid_matches:
        return False, "No tetradecanoic acid fragment found"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Check if at least one ester group is between tetradecanoic acid and alcohol/phenol
    for match in ester_matches:
        ester_oxygen_idx = match[0]
        ester_carbonyl_idx = match[1]
        ester_alcohol_idx = [atom.GetIdx() for atom in mol.GetAtomWithIdx(ester_carbonyl_idx).GetNeighbors()
                             if atom.GetIdx() != ester_oxygen_idx][0]

        # Check if the alcohol/phenol is connected to tetradecanoic acid
        if any(mol.GetAtomWithIdx(ester_alcohol_idx).GetNeighbors()[0].GetIdx() in tetradecanoic_acid_matches):
            return True, "Ester bond between tetradecanoic acid and alcohol/phenol found"

    return False, "No ester bond between tetradecanoic acid and alcohol/phenol found"