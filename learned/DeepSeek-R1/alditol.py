"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Get all carbon atoms
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    # Minimum 3 carbons (n ≥ 1 in general formula)
    if len(carbons) < 3:
        return False, "At least three carbons required"

    # Check linear chain structure
    end_carbons = 0
    for atom in carbons:
        # Count carbon neighbors for branching check
        carbon_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
        
        if carbon_neighbors == 1:
            end_carbons += 1
        elif carbon_neighbors != 2:
            return False, "Branched or non-linear chain"

    # Must have exactly two end carbons (linear chain)
    if end_carbons != 2:
        return False, "Not a linear chain"

    # Check hydroxyl groups
    for atom in carbons:
        hydroxyl_count = 0
        for neighbor in atom.GetNeighbors():
            # Check for -OH groups (oxygen with single bond and one hydrogen)
            if neighbor.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                    if neighbor.GetTotalNumHs() == 1:
                        hydroxyl_count += 1
        if hydroxyl_count != 1:
            return False, f"Carbon {atom.GetIdx()+1} has {hydroxyl_count} hydroxyl groups"

    return True, "Linear unbranched polyol with one hydroxyl per carbon"