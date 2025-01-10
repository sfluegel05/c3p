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

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found"
    
    # Check for small alcohol esters (methyl/ethyl) which are not waxes
    small_ester_pattern = Chem.MolFromSmarts("[CH3,CH2][OX2][CX3](=[OX1])")
    if mol.HasSubstructMatch(small_ester_pattern):
        # Verify if there's only a small ester (not part of a larger structure)
        if len(mol.GetSubstructMatches(small_ester_pattern)) == len(ester_matches):
            return False, "Simple methyl/ethyl ester, not a wax"

    # Count carbons and check for minimum chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:  # Typical waxes have at least 20 carbons total
        return False, f"Too few carbons ({c_count}) for a wax"

    # Look for long carbon chains (at least 8 carbons in a row)
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(long_chain_matches) < 1:
        return False, "Missing long carbon chains"

    # Count rotatable bonds to verify flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Not enough rotatable bonds for a wax"

    # Check molecular weight - waxes typically >350 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, "Molecular weight too low for a wax"

    # Check for carboxylic acids without ester groups
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(acid_pattern):
        return False, "Contains free carboxylic acid group"

    # Additional check for excessive branching - waxes are typically linear
    branching_pattern = Chem.MolFromSmarts("[CH1X4,CH0X4]")  # Carbons with 3 or 4 carbon neighbors
    branching_matches = mol.GetSubstructMatches(branching_pattern)
    if len(branching_matches) > 3:  # Allow some branching but not too much
        return False, "Too much branching for typical wax structure"

    # Check for long chain on both sides of ester
    for match in ester_matches:
        # Get the atoms connected to the ester oxygen
        ester_o = mol.GetAtomWithIdx(match[0])
        neighbors = [n.GetIdx() for n in ester_o.GetNeighbors()]
        for n in neighbors:
            if n != match[1]:  # not the carbonyl carbon
                # Check if this is part of a long chain
                if not any(n in chain for chain in long_chain_matches):
                    return False, "Ester group not connected to long chains on both sides"

    return True, "Long-chain ester with appropriate molecular weight and chain length"