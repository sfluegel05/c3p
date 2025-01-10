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

    # More flexible indole core patterns
    indole_patterns = [
        Chem.MolFromSmarts("[nH]1ccc2c1cccc2"),  # Basic indole
        Chem.MolFromSmarts("[nH]1ccc2c1cc([*])cc2"),  # Substituted indole
        Chem.MolFromSmarts("[nH]1ccc2c1c([*])ccc2"),  # Substituted indole
        Chem.MolFromSmarts("[nH]1ccc2c1c3ccccc32"),  # Fused ring system
        Chem.MolFromSmarts("[nH]1ccc2c1c3cc([*])cc32"),  # Fused and substituted
        Chem.MolFromSmarts("[nH]1ccc2c1c3c([*])cccc32"),  # Fused and substituted
        Chem.MolFromSmarts("[nH]1ccc2c1c3c([*])cc([*])c32")  # Fused and multiple substitutions
    ]
    
    has_indole = any(mol.HasSubstructMatch(pattern) for pattern in indole_patterns)
    if not has_indole:
        return False, "No indole core found"

    # More sophisticated terpenoid feature detection
    # Look for branched or cyclic carbon structures with at least 8 carbons
    terpenoid_patterns = [
        Chem.MolFromSmarts("[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]"),  # Linear
        Chem.MolFromSmarts("[C;!$(C=O)]1~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]1"),  # Cyclic
        Chem.MolFromSmarts("[C;!$(C=O)](~[C;!$(C=O)])(~[C;!$(C=O)])~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]"),  # Branched
        Chem.MolFromSmarts("[C;!$(C=O)]1~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]1"),  # Larger cyclic
        Chem.MolFromSmarts("[C;!$(C=O)](~[C;!$(C=O)])(~[C;!$(C=O)])~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]~[C;!$(C=O)]")  # Larger branched
    ]
    
    has_terpenoid = any(mol.HasSubstructMatch(pattern) for pattern in terpenoid_patterns)
    if not has_terpenoid:
        return False, "Insufficient terpenoid-like features"

    # Check for nitrogen atoms (from tryptophan precursor)
    nitrogen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if nitrogen_count < 1:
        return False, "Insufficient nitrogen atoms"

    # Check molecular weight range (typical for this class)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1200:
        return False, "Molecular weight out of expected range"

    # Check for ester groups (common in monoterpenoid indole alkaloids)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_count = len(mol.GetSubstructMatches(ester_pattern))
    if ester_count < 1:
        return False, "No ester groups found"

    # Check for connectivity between indole and terpenoid parts
    # Allow for more flexible connections through various functional groups
    linker_patterns = [
        Chem.MolFromSmarts("[nH]1ccc2c1cccc2~[*]~[C;!$(C=O)]"),  # Direct connection
        Chem.MolFromSmarts("[nH]1ccc2c1cccc2~[*]~[O,N]~[*]~[C;!$(C=O)]"),  # Through O or N
        Chem.MolFromSmarts("[nH]1ccc2c1cccc2~[*]~[C;!$(C=O)]~[*]~[C;!$(C=O)]"),  # Through C chain
        Chem.MolFromSmarts("[nH]1ccc2c1cccc2~[*]~[O,N]~[*]~[C;!$(C=O)]~[*]~[C;!$(C=O)]"),  # Longer chain
        Chem.MolFromSmarts("[nH]1ccc2c1cccc2~[*]~[C;!$(C=O)]~[*]~[O,N]~[*]~[C;!$(C=O)]")  # Alternate connection
    ]
    
    has_linker = any(mol.HasSubstructMatch(pattern) for pattern in linker_patterns)
    if not has_linker:
        return False, "No clear connection between indole and terpenoid parts"

    return True, "Contains indole core with terpenoid-like features, ester groups, and appropriate connectivity"