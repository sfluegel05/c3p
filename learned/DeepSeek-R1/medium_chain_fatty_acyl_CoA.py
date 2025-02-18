"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:61904 medium-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA has a CoA moiety linked via a thioester bond to a fatty acid with 6-12 carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for thioester group (C(=O)S-C) connected to CoA's sulfur
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2][CX4]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"

    # Verify CoA structure presence (adenine with amino group and multiple phosphates)
    # Adenine pattern matches 6-aminopurine core
    adenine_pattern = Chem.MolFromSmarts("n1c(N)nc2c1ncn2")
    # Phosphate pattern matches any P with double bond to O and at least two O neighbors
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)([OX2])[OX2]")
    
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine component"
    
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need at least 2"

    # Identify the fatty acid chain attached to the thioester
    try:
        # Get carbonyl carbon (attached to sulfur via thioester)
        sulfur_idx = thioester_matches[0][2]
        carbonyl_carbon = [a.GetIdx() for a in mol.GetAtomWithIdx(sulfur_idx).GetNeighbors() 
                         if a.GetSymbol() == "C" and a.GetBondToAtomIdx(sulfur_idx).GetBondType() == Chem.BondType.SINGLE][0]
    except IndexError:
        return False, "Invalid thioester connectivity"

    # Traverse the carbon chain from the carbonyl carbon (excluding CoA part)
    chain = []
    visited = set()
    
    def traverse(atom_idx, prev_idx):
        if atom_idx in visited:
            return
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        # Follow carbons and oxygen (for possible branching with ester groups)
        if atom.GetSymbol() in ["C", "O"]:
            if atom.GetSymbol() == "C":
                chain.append(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() != prev_idx:
                    traverse(neighbor.GetIdx(), atom_idx)
    
    traverse(carbonyl_carbon, sulfur_idx)

    # Calculate chain length (including carbonyl carbon)
    chain_length = len([a for a in chain if mol.GetAtomWithIdx(a).GetSymbol() == "C"])
    
    if not (6 <= chain_length <= 12):
        return False, f"Fatty acid chain length {chain_length} not in 6-12 range"

    return True, f"Medium-chain ({chain_length} carbons) fatty acyl-CoA with valid CoA structure"