"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: CHEBI:27774 chlorophyll
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from typing import Tuple

def is_chlorophyll(smiles: str) -> Tuple[bool, str]:
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophylls are a family of magnesium porphyrins with a fifth ring beyond the four pyrrole-like rings.
    The rings can have various side chains, typically including a long phytol chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for magnesium atom
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if len(mg_atoms) != 1:
        return False, "Chlorophylls must have exactly one magnesium atom"

    # Check for porphyrin ring system
    porphyrin_pattern = Chem.MolFromSmarts("[N;R]1[C;R]2[N;R]3[C;R]4[N;R]5[C;R]6[N;R]7[C;R]8[C;R]9[C;R]%10[C;R]%11[C;R]%12[C;R]%13[C;R]%14[C;R]%15[C;R]%16[C;R]%17[C;R]%18[C;R]%19[C;R]%20[C;R]1%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40")
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin ring system found"

    # Check for fifth ring beyond the four pyrrole-like rings
    fifth_ring_pattern = Chem.MolFromSmarts("[N;R]1[C;R]2[N;R]3[C;R]4[N;R]5[C;R]6[N;R]7[C;R]8[C;R]9[C;R]%10[C;R]%11[C;R]%12[C;R]%13[C;R]%14[C;R]%15[C;R]%16[C;R]%17[C;R]%18[C;R]%19[C;R]%20[C;R]1%21%22%23%24%25%26%27%28%29%30%31%32%33%34%35%36%37%38%39%40[C;R]%41[C;R]%42[C;R]%43")
    if not mol.HasSubstructMatch(fifth_ring_pattern):
        return False, "No fifth ring found beyond the four pyrrole-like rings"

    # Check for long aliphatic side chain (phytol or similar)
    phytol_pattern = Chem.MolFromSmarts("[CH3][CX4]([CH3])[CH2][CH2][CH2][CX4]([CH3])[CH2][CH2][CH2][CH2][CX4]([CH3])[CH2][CH2][CH2][CX4]([CH3])[CH2][CH2][CH2][CX4]([CH3])[CH2][CH2][CH2][CX4]([CH3])[CH2][CH2][CH2][CH2][CH2][CH2]")
    phytol_matches = mol.GetSubstructMatches(phytol_pattern)
    if not phytol_matches:
        return False, "No long aliphatic side chain (phytol or similar) found"

    # Additional checks for typical chlorophyll substituents
    vinyl_pattern = Chem.MolFromSmarts("[CH2]=[CH]")
    vinyl_matches = mol.GetSubstructMatches(vinyl_pattern)

    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    ketone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[#6]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    if not (vinyl_matches and ester_matches and ketone_matches):
        return False, "Missing typical chlorophyll substituents (vinyl, ester, ketone)"

    return True, "Contains a porphyrin ring system with a magnesium atom, a fifth ring, and a long aliphatic side chain (phytol or similar). Typical chlorophyll substituents (vinyl, ester, ketone) are also present."