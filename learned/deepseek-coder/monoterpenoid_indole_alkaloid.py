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

    # More flexible indole core pattern that allows substitutions and fused rings
    indole_patterns = [
        Chem.MolFromSmarts("[nH]1ccc2c1cccc2"),  # Basic indole
        Chem.MolFromSmarts("[nH]1ccc2c1cc([*])cc2"),  # Substituted indole
        Chem.MolFromSmarts("[nH]1ccc2c1c([*])ccc2"),  # Substituted indole
        Chem.MolFromSmarts("[nH]1ccc2c1c3ccccc32")  # Fused ring system
    ]
    
    has_indole = any(mol.HasSubstructMatch(pattern) for pattern in indole_patterns)
    if not has_indole:
        return False, "No indole core found"

    # More sophisticated terpenoid feature detection
    # Look for at least 8 carbons in a branched or cyclic arrangement
    terpenoid_pattern = Chem.MolFromSmarts("[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]")
    if not mol.HasSubstructMatch(terpenoid_pattern):
        return False, "Insufficient terpenoid-like features"

    # Check for nitrogen atoms (from tryptophan precursor)
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 1:
        return False, "Insufficient nitrogen atoms"

    # Check molecular weight range (typical for this class)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1200:
        return False, "Molecular weight out of expected range"

    # Check for connectivity between indole and terpenoid parts
    # Look for a carbon atom connected to both systems
    linker_pattern = Chem.MolFromSmarts("[nH]1ccc2c1cccc2~[*]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]")
    if not mol.HasSubstructMatch(linker_pattern):
        return False, "No clear connection between indole and terpenoid parts"

    # Check for appropriate ring count (typically 2-6 rings)
    ring_count = len(mol.GetRingInfo().AtomRings())
    if ring_count < 2 or ring_count > 6:
        return False, "Ring count out of expected range"

    return True, "Contains indole core with terpenoid-like features and appropriate connectivity"