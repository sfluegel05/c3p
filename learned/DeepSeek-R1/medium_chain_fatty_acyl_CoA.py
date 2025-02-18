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

    # Check for thioester group (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"

    # Verify CoA structure: adenine, phosphopantetheine, and phosphates
    # Adenine pattern (6-aminopurine)
    adenine_pattern = Chem.MolFromSmarts("n1c(nc2c1ncn2)N")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Missing adenine component"

    # CoA typically has 3 phosphate groups (diphosphate bridge + phosphopantetheine)
    phosphate_pattern = Chem.MolFromSmarts("[PX4](=O)([OX2])[OX2]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 3:
        return False, f"Found {len(phosphate_matches)} phosphate groups, need at least 3"

    # Identify the fatty acid chain attached to the thioester
    try:
        # Get sulfur atom from thioester
        sulfur_idx = thioester_matches[0][2]
        # Get carbonyl carbon (attached to sulfur)
        carbonyl_carbon = [n.GetIdx() for n in mol.GetAtomWithIdx(sulfur_idx).GetNeighbors() 
                          if n.GetSymbol() == "C" and 
                          mol.GetBondBetweenAtoms(n.GetIdx(), sulfur_idx).GetBondType() == Chem.BondType.SINGLE][0]
    except (IndexError, AttributeError):
        return False, "Invalid thioester connectivity"

    # Traverse the carbon chain from carbonyl carbon (excluding CoA part)
    chain = []
    visited = set()
    
    def traverse(atom_idx, prev_idx):
        if atom_idx in visited:
            return
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == "C":
            chain.append(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() != prev_idx:
                    traverse(neighbor.GetIdx(), atom_idx)
    
    traverse(carbonyl_carbon, sulfur_idx)

    # Calculate chain length (number of carbons in fatty acid R group)
    # Add 1 to include the carbonyl carbon in total fatty acid length
    chain_length = len(chain)
    total_fa_length = chain_length + 1  # R group carbons + carbonyl
    
    if not (6 <= total_fa_length <= 12):
        return False, f"Fatty acid chain length {total_fa_length} not in 6-12 range"

    return True, f"Medium-chain ({total_fa_length} carbons) fatty acyl-CoA with valid CoA structure"