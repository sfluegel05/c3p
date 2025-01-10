"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:76224 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.
    A short-chain fatty acyl-CoA is a fatty acyl-CoA with a short-chain fatty acid (2-6 carbons)
    attached to CoA via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety using a more flexible pattern
    # This pattern matches the core structure of CoA, allowing for different protonation states
    coa_pattern = Chem.MolFromSmarts("[O-]P(=O)([O-])OP(=O)([O-])OCC1OC([C@H](O)[C@@H]1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Check for thioester bond (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) == 0:
        return False, "No thioester bond found"

    # Extract the fatty acid chain attached to the thioester
    # The fatty acid chain should be 2-6 carbons long
    fatty_acid_chain = False
    for match in thioester_matches:
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == "C" and atom.GetDegree() == 3:  # Carbon in thioester bond
                # Traverse the chain to count carbons
                chain_length = 0
                stack = [(atom, 0)]
                visited = set()
                while stack:
                    current_atom, depth = stack.pop()
                    if current_atom.GetIdx() in visited:
                        continue
                    visited.add(current_atom.GetIdx())
                    if current_atom.GetSymbol() == "C":
                        chain_length += 1
                    if chain_length > 6:
                        break
                    for neighbor in current_atom.GetNeighbors():
                        if neighbor.GetSymbol() != "S" and neighbor.GetSymbol() != "O":
                            stack.append((neighbor, depth + 1))
                if 2 <= chain_length <= 6:
                    fatty_acid_chain = True
                    break
        if fatty_acid_chain:
            break

    if not fatty_acid_chain:
        return False, "Fatty acid chain length not in the range of 2-6 carbons"

    return True, "Contains CoA moiety with a short-chain fatty acid (2-6 carbons) attached via a thioester bond"