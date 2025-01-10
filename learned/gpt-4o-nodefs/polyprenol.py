"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    Polyprenols are oligomers with three or more isoprene units, often ending in
    an alcohol group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for isoprene unit: C=C-C-C=C
    isoprene_smarts = "C=C-C-C=C"
    isoprene_pattern = Chem.MolFromSmarts(isoprene_smarts)

    # Check for repeated isoprene units - at least three times (polyprenol definition)
    if mol.HasSubstructMatch(isoprene_pattern):
        if len(mol.GetSubstructMatches(isoprene_pattern)) >= 3:
            # Verify end group contains alcohol (OH)
            if mol.GetAtoms()[-1].GetSymbol() == "O" and mol.GetAtoms()[-2].GetSymbol() == "C":
                return True, "Molecule is a polyprenol with repeating isoprene units ending in an alcohol group"
            else:
                return False, "Isoprene units are present, but molecule does not end in an alcohol group"
        else:
            return False, f"Only {len(mol.GetSubstructMatches(isoprene_pattern))} isoprene units found, at least 3 required"
    else:
        return False, "No isoprene units found"