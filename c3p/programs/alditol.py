"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:17522 alditol
"""
from rdkit import Chem
from rdkit.Chem import Mol

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

    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings"

    # Check for branching (each carbon has <=2 carbon neighbors)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            carbon_neighbors = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetAtomicNum() == 6)
            if carbon_neighbors > 2:
                return False, "Branched carbon chain"

    # Check all carbons have at least one hydroxyl group
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            has_hydroxyl = any(neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1 
                              for neighbor in atom.GetNeighbors())
            if not has_hydroxyl:
                return False, f"Carbon {atom.GetIdx()+1} has no hydroxyl group"

    return True, "Acyclic, unbranched polyol with hydroxyls on all carbons"