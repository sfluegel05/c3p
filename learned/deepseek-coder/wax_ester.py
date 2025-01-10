"""
Classifies: CHEBI:10036 wax ester
"""
"""
Classifies: CHEBI:17855 wax ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of a fatty acid with a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group (-COO-)
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for long hydrocarbon chains (at least 8 carbons in each chain)
    # We look for a pattern of at least 8 consecutive carbons, allowing for double bonds and branching
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(long_chain_matches) < 2:
        return False, f"Found {len(long_chain_matches)} long chains, need at least 2"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Chains too short to be fatty acid and fatty alcohol"

    # Check molecular weight - wax esters typically >200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for wax ester"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 16:
        return False, "Too few carbons for wax ester"
    if o_count < 2:
        return False, "Must have at least 2 oxygens (1 ester group)"

    # Exclude molecules with additional functional groups or structures
    # Exclude cholesteryl esters, retinyl esters, etc.
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)OC(=O))C)C")):
        return False, "Excluded cholesteryl ester"
    if mol.HasSubstructMatch(Chem.MolFromSmarts("C1=C(C(=O)C=C(C1=O)C)C")):
        return False, "Excluded retinyl ester"
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4][OX2][CX3](=[OX1])")):
        return False, "Excluded molecules with multiple ester groups"

    return True, "Contains one ester group with two long hydrocarbon chains"