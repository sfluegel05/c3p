"""
Classifies: CHEBI:33447 phospho sugar
"""
"""
Classifies: CHEBI:27838 phospho sugar
A phospho sugar is any monosaccharide containing an alcoholic hydroxy group esterified with phosphoric acid.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for phosphate group patterns
    phosphate_pattern = Chem.MolFromSmarts("[P+](O)(O)(O)=O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Look for monosaccharide backbone patterns
    monosaccharide_patterns = [Chem.MolFromSmarts("[OX2]r1[CX4]([OX2])[CX4]([OX2])[CX4]([OX2])[CX4]([OX2])1"),  # Pyranose ring
                               Chem.MolFromSmarts("[OX2]r1[CX4]([OX2])[CX4]([OX2])[CX4]([OX2])1")]  # Furanose ring

    monosaccharide_match = False
    for pattern in monosaccharide_patterns:
        if len(mol.GetSubstructMatches(pattern)) == 1:
            monosaccharide_match = True
            break

    if not monosaccharide_match:
        return False, "No monosaccharide backbone found"

    # Check for phosphate-sugar connectivity
    phospho_sugar_pattern = Chem.MolFromSmarts("[P+](O)(O)(O)=O~[OX2]~[CX4]")
    if not mol.HasSubstructMatch(phospho_sugar_pattern):
        return False, "Phosphate group not connected to an alcoholic hydroxy group"

    # Additional checks
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 3 or c_count > 8:
        return False, "Number of carbons outside the expected range for monosaccharides"

    if o_count < 3 or o_count > 8:
        return False, "Number of oxygens outside the expected range for monosaccharides"

    return True, "Contains a phosphate group attached to an alcoholic hydroxy group of a monosaccharide"