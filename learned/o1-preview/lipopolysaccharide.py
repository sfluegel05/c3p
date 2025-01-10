"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: lipopolysaccharide
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    A lipopolysaccharide consists of saccharide units attached to long-chain fatty acids
    via ester or amide linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for sugar rings (pyranose and furanose forms)
    sugar_smarts = Chem.MolFromSmarts("[C;R][$([O][C;R])][C;R][$([O][C;R])][C;R][$([O][C;R])][C;R]")  # Simplified sugar pattern
    sugar_matches = mol.GetSubstructMatches(sugar_smarts)
    if len(sugar_matches) == 0:
        return False, "No sugar units found"

    # Look for long-chain fatty acids (aliphatic chains with >10 carbons)
    fatty_acid_smarts = Chem.MolFromSmarts("C(=O)[O,N][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4][C;X4]")  # Chain of at least 10 carbons
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_smarts)
    if len(fatty_acid_matches) == 0:
        return False, "No long-chain fatty acids found"

    # Check for ester or amide linkages between sugars and fatty acids
    ester_or_amide_linkage = Chem.MolFromSmarts("[C;R][$([O][C;R])][C;R][O,N][C](=O)[C;X4]")
    linkage_matches = mol.GetSubstructMatches(ester_or_amide_linkage)
    if len(linkage_matches) == 0:
        return False, "No ester or amide linkages between sugars and fatty acids found"

    # Check for overall size - lipopolysaccharides are large molecules
    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 1000:
        return False, "Molecular weight too low for lipopolysaccharide"

    return True, "Contains sugar units linked to long-chain fatty acids via ester or amide bonds"

__metadata__ = {
    'chemical_class': {
        'name': 'lipopolysaccharide',
        'definition': 'Liposaccharide natural compounds consisting of a trisaccharide repeating unit (two heptose units and octulosonic acid) with oligosaccharide side chains and 3-hydroxytetradecanoic acid units (they are a major constituent of the cell walls of Gram-negative bacteria).'
    }
}