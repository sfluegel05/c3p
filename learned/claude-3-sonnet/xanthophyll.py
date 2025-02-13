"""
Classifies: CHEBI:27325 xanthophyll
"""
"""
Classifies: CHEBI:38808 xanthophyll
A subclass of carotenoids consisting of the oxygenated carotenes.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_xanthophyll(smiles: str):
    """
    Determines if a molecule is a xanthophyll based on its SMILES string.
    A xanthophyll is a carotenoid containing oxygen atoms in the form of hydroxyl,
    epoxy, or keto groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a xanthophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carotenoid backbone pattern
    carotenoid_pattern = Chem.MolFromSmarts("[C;R]=[C;R][C;R]=[C;R][C;R]=[C;R][C;R]=[C;R][C;R]=[C;R]")
    if not mol.HasSubstructMatch(carotenoid_pattern):
        return False, "No carotenoid backbone found"

    # Look for oxygen atoms
    oxygen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    if not oxygen_atoms:
        return False, "No oxygen atoms found"

    # Check for specific functional groups and positions
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    epoxy_pattern = Chem.MolFromSmarts("[O;R]")
    keto_pattern = Chem.MolFromSmarts("[C=O]")

    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    epoxy_matches = mol.GetSubstructMatches(epoxy_pattern)
    keto_matches = mol.GetSubstructMatches(keto_pattern)

    # Xanthophylls typically have oxygen atoms at specific positions (3, 3', 4, 4', 5, 6, 5', 6')
    oxygen_positions = [mol.GetAtomWithIdx(match[0]).GetProp("nAtomLabel") for matches in [hydroxy_matches, epoxy_matches, keto_matches] for match in matches]
    if any(pos in ["3", "3'", "4", "4'", "5", "6", "5'", "6'"] for pos in oxygen_positions):
        return True, "Contains a carotenoid backbone with oxygen atoms (hydroxy, epoxy, or keto groups) at typical xanthophyll positions"
    else:
        return False, "Oxygen atoms not present at typical xanthophyll positions"