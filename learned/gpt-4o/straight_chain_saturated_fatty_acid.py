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
        bool: True if molecule is a straight-chain saturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (must be terminal)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    if len(carboxylic_matches) != 1:
        return False, "No terminal carboxylic acid group found or multiple groups detected"

    # Ensure carbon chain saturation (apart from carboxylic)
    # Saturation refers to carbon-carbon bonds only
    for bond in mol.GetBonds():
        if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return False, "Presence of non-single C-C bonds indicates unsaturation"

    # Determine if molecule is a straight chain
    # Count carbon atoms; ensure they are forming a continuous chain
    c_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(c_atoms) < 3:
        return False, "Too few carbon atoms for a fatty acid"

    # We require a linear alignment with carboxylic acid
    terminal_carbon_index = carboxylic_matches[0][0]
    carbon_chain = [terminal_carbon_index]
    current_index = terminal_carbon_index
    
    while True:
        # Extend the chain forward
        next_carbon = None
        for neighbor in mol.GetAtomWithIdx(current_index).GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in carbon_chain:
                next_carbon = neighbor.GetIdx()
                break
        
        if next_carbon is None:
            break
        
        carbon_chain.append(next_carbon)
        current_index = next_carbon

    if len(carbon_chain) != len(c_atoms):
        return False, "Not a continuous straight chain; branches or side chains detected"

    # Hydroxy group allowance: count OH groups apart from in carboxylic acid
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)

    if any(oxyg_match[0] in chain for chain in hydroxy_matches):
        if len(hydroxy_matches) > 1:
            return False, f"Too many hydroxy groups: found {len(hydroxy_matches)}, need at most 1"

    return True, "Molecule is a straight-chain saturated fatty acid"

# Example use:
# smiles = "CCCCCCCC(O)=O"  # Octanoic acid
# result, reason = is_straight_chain_saturated_fatty_acid(smiles)
# print(result, reason)