"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: CHEBI:35195 wax
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is typically a long-chain ester formed between a fatty acid and a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for exactly one ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1 for wax"

    # Count carbons and check for minimum chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:  # Typical waxes have at least 20 carbons total
        return False, f"Too few carbons ({c_count}) for a wax"

    # Look for long carbon chains
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 2:  # Need at least 2 long chains (one from alcohol, one from acid)
        return False, "Missing long carbon chains"

    # Count rotatable bonds to verify flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Not enough rotatable bonds for a wax"

    # Check molecular weight - waxes typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for a wax"

    # Count oxygens - waxes typically have 2 oxygens (from single ester group)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 2:
        return False, f"Found {o_count} oxygens, waxes typically have exactly 2 (ester group)"

    # Additional check for branching - waxes are typically linear
    branching_pattern = Chem.MolFromSmarts("[CH1X4,CH0X4]")  # Carbons with 3 or 4 carbon neighbors
    branching_matches = mol.GetSubstructMatches(branching_pattern)
    if len(branching_matches) > 2:  # Allow some branching but not too much
        return False, "Too much branching for typical wax structure"

    return True, "Long-chain ester with appropriate molecular weight and chain length"