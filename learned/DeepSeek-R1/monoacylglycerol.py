"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: CHEBI:17478 monoacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol has a glycerol backbone with exactly one acyl (ester) group
    and two other -OH or alkyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone (C-O-C-O-C pattern with three carbons)
    glycerol_pattern = Chem.MolFromSmarts("[CH2]-[CH](-[CH2])-O")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Find all ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check that the ester is attached to the glycerol backbone
    # Get the oxygen atom in the ester
    ester_o = ester_matches[0][0]
    ester_neighbors = [x.GetIdx() for x in mol.GetAtomWithIdx(ester_o).GetNeighbors() if x.GetSymbol() == 'C']
    if not ester_neighbors:
        return False, "Ester not attached to carbon chain"
    
    # Verify the ester is connected to a carbon chain (at least 2 carbons)
    chain_length = 0
    visited = set()
    stack = [(ester_neighbors[0], 0)]
    while stack:
        atom_idx, depth = stack.pop()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() == 'C':
            chain_length = max(chain_length, depth + 1)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() not in visited:
                    stack.append((neighbor.GetIdx(), depth + 1))
    if chain_length < 2:
        return False, "Acyl chain too short"

    # Count hydroxyl groups (should be two if others are H, or check for ethers)
    hydroxyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            if atom.GetTotalNumHs() >= 1 and atom.GetDegree() == 1:
                hydroxyl_count += 1
    if hydroxyl_count < 2:
        # Check for ethers (O connected to two carbons)
        ether_o = 0
        for atom in mol.GetAtoms():
            if atom.GetSymbol() == 'O' and atom.GetDegree() == 2:
                ether_o += 1
        if ether_o + hydroxyl_count < 2:
            return False, "Insufficient hydroxyl/ether groups on glycerol"

    return True, "Contains glycerol backbone with one acyl group and two hydroxyl/alkyl groups"