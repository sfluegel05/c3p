"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the carboxylic acid group pattern
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylic acid group found"

    # Make sure there are no double or triple bonds in the carbon chain
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False, "Carbon chain contains double or triple bonds"
    
    # Check if there are branches - start from carboxyl group
    chain_lengths = []
    for match in mol.GetSubstructMatches(carboxylate_pattern):
        carboxyl_carbon = mol.GetAtomWithIdx(match[0])
        visited = set()
        stack = [(carboxyl_carbon, None, 0)]  # (atom, previous atom, chain length)
        
        while stack:
            atom, prev_atom, length = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())

            # If it's a terminal carbon, record the chain length
            if atom.GetAtomicNum() == 6 and len([b for b in atom.GetBonds() if b.GetOtherAtom(atom).GetAtomicNum() == 6]) == 1:
                chain_lengths.append(length + 1)

            # Add neighboring carbons to check
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetIdx() != (prev_atom.GetIdx() if prev_atom else None) and neighbor.GetAtomicNum() == 6:
                    stack.append((neighbor, atom, length+1))

    if len(chain_lengths) > 1:
        return False, "Carbon chain is branched"

    if not chain_lengths:
        return False, "No valid carbon chain found"

    return True, "Molecule is a straight-chain saturated fatty acid"

# Example for testing
example_smiles = "CCCCCCC(O)=O"  # Heptanoic acid
result, reason = is_straight_chain_saturated_fatty_acid(example_smiles)
print(f"{example_smiles}: {result} ({reason})")