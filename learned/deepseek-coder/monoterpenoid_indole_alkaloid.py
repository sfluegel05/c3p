"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: CHEBI:71934 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    These compounds are biosynthesized from L-tryptophan and diisoprenoid building blocks.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for indole core pattern
    indole_pattern = Chem.MolFromSmarts("[nH]1ccc2c1cccc2")
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole core found"

    # Check for terpenoid-like features (isoprene units)
    # Look for at least 2 isoprene-like patterns (C5H8 units)
    isoprene_pattern = Chem.MolFromSmarts("CC(=C)C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, "Insufficient terpenoid-like features"

    # Check for nitrogen atoms (from tryptophan precursor)
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 2:
        return False, "Insufficient nitrogen atoms"

    # Check molecular weight range (typical for this class)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1000:
        return False, "Molecular weight out of expected range"

    # Check for connectivity between indole and terpenoid parts
    # Look for a carbon atom connected to both systems
    linker_pattern = Chem.MolFromSmarts("[nH]1ccc2c1cccc2.[*]~C~[*]")
    if not mol.HasSubstructMatch(linker_pattern):
        return False, "No clear connection between indole and terpenoid parts"

    # Check for typical functional groups (esters, ethers, etc.)
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    ether_pattern = Chem.MolFromSmarts("[#6][OX2][#6]")
    if not (mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(ether_pattern)):
        return False, "Missing typical functional groups"

    # Check for appropriate ring count (typically 3-5 rings)
    ring_count = len(mol.GetRingInfo().AtomRings())
    if ring_count < 3 or ring_count > 5:
        return False, "Ring count out of expected range"

    return True, "Contains indole core with terpenoid-like features and appropriate connectivity"