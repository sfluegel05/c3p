"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: Aliphatic Alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is an alcohol derived from an aliphatic compound.
    This function checks if the molecule has at least one hydroxyl group (-OH)
    attached to a non-aromatic, sp3-hybridized carbon atom (aliphatic alcohol group).
    It also excludes molecules that are likely sugars or peptides.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for hydroxyl group attached to non-aromatic, sp3-hybridized carbon
    alcohol_pattern = Chem.MolFromSmarts('[CX4;!$(C[-])(H)][OX2H]')

    # Find matches of the alcohol pattern
    matches = mol.GetSubstructMatches(alcohol_pattern)
    if not matches:
        return False, "No aliphatic alcohol groups found"

    # Exclude molecules that are likely sugars (glycosides)
    # Sugars often have multiple hydroxyl groups attached to sp3 carbons and ring structures
    # Check for glycosidic bonds (C-O-C linkages in rings)
    glycosidic_pattern = Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H]1O')
    if mol.HasSubstructMatch(glycosidic_pattern):
        return False, "Molecule is likely a sugar or glycoside"

    # Exclude molecules that are peptides (contain peptide bonds)
    # Peptide bond pattern: C(=O)N
    peptide_pattern = Chem.MolFromSmarts('C(=O)N')
    if mol.HasSubstructMatch(peptide_pattern):
        return False, "Molecule is likely a peptide"

    # If molecule passes the above checks, it is considered an aliphatic alcohol
    return True, "Contains aliphatic alcohol group(s)"