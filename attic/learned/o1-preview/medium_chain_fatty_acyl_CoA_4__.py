"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    A medium-chain fatty acyl-CoA(4-) is an acyl-CoA molecule with a fatty acyl chain
    of 6 to 12 carbons attached via a thioester bond to Coenzyme A, and with deprotonated
    phosphate groups resulting in a 4- charge at physiological pH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for adenine moiety
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2N")
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, "Adenine moiety not found"

    # Look for thioester bond (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "Thioester bond not found"

    # Get the acyl chain length
    acyl_carbon_idx = thioester_matches[0][0]  # Carbonyl carbon of thioester
    sulfur_idx = thioester_matches[0][2]       # Sulfur atom of thioester bond

    # Get the acyl chain starting atom (carbon attached to carbonyl carbon excluding O and S)
    carbonyl_carbon = mol.GetAtomWithIdx(acyl_carbon_idx)
    acyl_chain_atom = None
    for neighbor in carbonyl_carbon.GetNeighbors():
        idx = neighbor.GetIdx()
        if idx != sulfur_idx and neighbor.GetAtomicNum() != 8:
            acyl_chain_atom = neighbor
            break
    if acyl_chain_atom is None:
        return False, "Could not find acyl chain attached to carbonyl carbon"

    # Traverse the acyl chain and count carbons
    visited = set()
    stack = [acyl_chain_atom]
    acyl_chain_length = 0

    while stack:
        atom = stack.pop()
        idx = atom.GetIdx()
        if idx in visited:
            continue
        visited.add(idx)
        if atom.GetAtomicNum() == 6:
            acyl_chain_length += 1
            for neighbor in atom.GetNeighbors():
                n_idx = neighbor.GetIdx()
                if neighbor.GetAtomicNum() == 6 and n_idx not in visited and neighbor.GetIdx() != carbonyl_carbon.GetIdx():
                    stack.append(neighbor)

    if acyl_chain_length < 6 or acyl_chain_length > 12:
        return False, f"Acyl chain length is {acyl_chain_length}, not in the range 6-12 for medium-chain fatty acids"

    # Count phosphorus atoms (phosphate groups)
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count < 3:
        return False, f"Found {p_count} phosphorus atoms, expected at least 3 phosphate groups in CoA"

    # Return success
    return True, "Molecule is a medium-chain fatty acyl-CoA(4-) with appropriate acyl chain length and CoA structure"