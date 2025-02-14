"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import rdMolStandardize

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    A glycosaminoglycan is a polysaccharide containing a substantial proportion of aminomonosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use SugarRemover to identify sugar units
    remover = rdMolStandardize.SugarRemover()
    sugar_fragments = remover.GetSugarFragments(mol)

    total_sugars = len(sugar_fragments)
    if total_sugars == 0:
        return False, "No sugar units found"

    # Identify amino sugars
    amino_sugar_count = 0
    for sugar in sugar_fragments:
        has_nitrogen = False
        for atom in sugar.GetAtoms():
            if atom.GetSymbol() == 'N':
                has_nitrogen = True
                break
        if has_nitrogen:
            amino_sugar_count += 1

    if amino_sugar_count == 0:
        return False, "No amino sugars found"

    # Check proportion of amino sugars
    proportion = amino_sugar_count / total_sugars
    if proportion < 0.3:
        return False, f"Only {amino_sugar_count} out of {total_sugars} sugars are amino sugars"

    # Check for polysaccharide chain (more than 3 sugar units)
    if total_sugars < 3:
        return False, f"Only {total_sugars} sugar units found, not enough for a polysaccharide"

    return True, f"Contains polysaccharide chain with {amino_sugar_count} amino sugars out of {total_sugars} sugar units"