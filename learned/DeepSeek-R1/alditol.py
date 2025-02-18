"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts

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
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings"
    
    # Check for carboxylic acid groups
    carboxylic_acid = MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(carboxylic_acid):
        return False, "Contains carboxylic acid group"
    
    # Get all carbon atoms
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) < 3:
        return False, "At least three carbons required"
    
    # Check linear chain structure (exactly two end carbons)
    end_carbons = 0
    for atom in carbons:
        carbon_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6)
        if carbon_neighbors == 1:
            end_carbons += 1
        elif carbon_neighbors != 2:
            return False, "Branched or non-linear chain"
    
    if end_carbons != 2:
        return False, "Not a linear chain"
    
    # Check each carbon has at least one hydroxyl and no more than one
    for atom in carbons:
        hydroxyl_count = 0
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                    # Check if it's a hydroxyl (either -OH or substituted -O- group)
                    # Count as hydroxyl if the oxygen has at least one H (could be part of ether)
                    # This is a simplification and may not be accurate
                    if neighbor.GetTotalNumHs() >= 1:
                        hydroxyl_count += 1
        if hydroxyl_count < 1:
            return False, f"Carbon {atom.GetIdx()+1} has no hydroxyl group"
        if hydroxyl_count > 1:
            return False, f"Carbon {atom.GetIdx()+1} has multiple hydroxyl groups"
    
    # Check ends are CH2OH (at least two Hs and one hydroxyl)
    for atom in carbons:
        if sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 6) == 1:  # end carbon
            h_count = atom.GetTotalNumHs()
            if h_count < 2:
                return False, f"End carbon {atom.GetIdx()+1} has insufficient hydrogens"
            hydroxyl_found = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond and bond.GetBondType() == Chem.BondType.SINGLE:
                        hydroxyl_found = True
            if not hydroxyl_found:
                return False, f"End carbon {atom.GetIdx()+1} lacks hydroxyl group"
    
    return True, "Linear acyclic polyol with one hydroxyl per carbon in main chain"