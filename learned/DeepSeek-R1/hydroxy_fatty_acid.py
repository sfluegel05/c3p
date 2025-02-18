"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: CHEBI hydroxy fatty acid
"""
from rdkit import Chem

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is a fatty acid (carboxylic acid with aliphatic chain)
    carrying one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find carboxylic acid groups (including deprotonated)
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1,O-]")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carboxylic_matches:
        return False, "No carboxylic acid group"
    carboxylic_carbons = {match[0] for match in carboxylic_matches}

    # Find hydroxyl groups not part of carboxylic acid
    hydroxyl_oxygens = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1:
            # Check if oxygen is in a hydroxyl group (single bond to carbon)
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 1 and neighbors[0].GetAtomicNum() == 6:
                # Verify not part of carboxylic acid group
                if not any(atom.GetIdx() in match for match in carboxylic_matches):
                    hydroxyl_oxygens.append(atom)

    if not hydroxyl_oxygens:
        return False, "No hydroxyl groups outside carboxylic acid"

    # Check if any hydroxyl is connected via aliphatic chain to carboxylic acid
    for o_atom in hydroxyl_oxygens:
        hydroxyl_carbon = o_atom.GetNeighbors()[0]
        for cc_idx in carboxylic_carbons:
            try:
                path = Chem.GetShortestPath(mol, hydroxyl_carbon.GetIdx(), cc_idx)
            except:
                continue  # No path exists
            
            # Verify all atoms in path are aliphatic carbons
            valid = True
            for atom_idx in path:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() != 6 or atom.GetIsAromatic():
                    valid = False
                    break
                # Check bonds in path (optional but adds robustness)
                if atom_idx != path[0]:
                    bond = mol.GetBondBetweenAtoms(prev_idx, atom_idx)
                    if bond.GetBondType() == Chem.BondType.AROMATIC:
                        valid = False
                        break
                prev_idx = atom_idx
            
            if valid:
                return True, "Hydroxyl group on aliphatic chain connected to carboxylic acid"

    return False, "No hydroxyl group on aliphatic chain connected to carboxylic acid"