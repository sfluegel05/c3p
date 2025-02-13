"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions as Reactions

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is a molecule with a glycerol moiety where the acyl group
    is esterified at position 1.

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

    # Glycerol backbone pattern: C(O)C(O)CO
    glycerol_pattern = Chem.MolFromSmarts("C(O)C(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Ester linkage pattern at primary alcohol position: OC(=O)C 
    ester_at_1_position = Chem.MolFromSmarts("COC(=O)C")
    if not mol.HasSubstructMatch(ester_at_1_position):
        return False, "No ester linkage at 1-position found"

    # Ensure there is only one acyl chain esterified (checks for only one ester linkage)
    ester_groups = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O"))
    if len(ester_groups) != 1:
        return False, f"Found {len(ester_groups)} ester groups, need exactly 1"

    # Glycerol specific chirality check (if Stereochemistry matters)
    chiral_center = Chem.MolFromSmiles("[C@@H]")
    if mol.HasSubstructMatch(chiral_center):
        chiral_matches = mol.GetSubstructMatches(chiral_center)
        if len(chiral_matches) == 0:
            return False, "Needs specific glycerol chirality"

    return True, "Contains glycerol backbone with single acyl group esterified at 1-position"