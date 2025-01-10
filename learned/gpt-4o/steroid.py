"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid usually has a cyclopenta[a]phenanthrene carbon skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Cyclopenta[a]phenanthrene skeleton pattern: 4 rings (3x6-membered, 1x5-membered)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C3C(C4=CC=CC=C4)C2C1C3")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not contain cyclopenta[a]phenanthrene skeleton"
    
    # Look for methyl groups, typically at C-10 and C-13
    methyl_pattern = Chem.MolFromSmarts("C(C)(C)C")
    methyl_matches = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' and len(atom.GetNeighbors()) == 4]
    if len(methyl_matches) < 2:
        return False, "Missing methyl groups at key positions"

    # Optional: Look for an alkyl group at C-17
    alkyl_pattern = Chem.MolFromSmarts("CC")
    if not mol.HasSubstructMatch(alkyl_pattern):
        return False, "Missing alkyl group"

    # Check stereochemistry: at key chiral centers
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if len(chiral_centers) < 4:
        return False, "Incorrect number of chiral centers"

    return True, "Contains features typical of a steroid: cyclopenta[a]phenanthrene skeleton with methyl groups"