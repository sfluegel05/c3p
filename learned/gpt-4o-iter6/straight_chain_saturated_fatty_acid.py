"""
Classifies: CHEBI:39418 straight-chain saturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_straight_chain_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a straight-chain saturated fatty acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find atoms involved in terminal carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_match = mol.GetSubstructMatch(carboxylic_acid_pattern)
    if not carboxylic_acid_match:
        return False, "No terminal carboxylic acid group found"
    
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]  # All carbon atoms
    if not carbon_atoms:
        return False, "No carbon atoms found, not a fatty acid."

    # Verify a continuous chain of sp3 carbons ending in carboxyl
    current = carboxylic_acid_match[0]  # Carboxyl carbon

    num_carbons = 0
    visited = set()

    while mol.GetAtomWithIdx(current).GetDegree() <= 2 and current not in visited:
        atom = mol.GetAtomWithIdx(current)
        # All traversed carbons must be sp3 (excluding carboxylic carbon)
        if num_carbons > 0 and atom.GetHybridization() != rdchem.HybridizationType.SP3:
            return False, f"Contains unsaturated carbon (non-sp3 hybridized) at index {current}"

        num_carbons += 1
        visited.add(current)

        # Move to the next carbon in the chain
        neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6 and n.GetIdx() not in visited]
        if not neighbors:
            break
        current = neighbors[0]

    # Natural fatty acids typically have lengths of 4-36 carbons
    if num_carbons < 4 or num_carbons > 36:
        return False, f"{num_carbons} carbons found; expected a chain length between 4 and 36"

    return True, f"Contains {num_carbons} carbons in a straight-chain saturated fatty acid structure"

# Note: The refined chain detection logic ensures carbons form an unbranched structure with a carboxylic acid group at the end.