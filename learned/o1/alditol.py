"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:37277 alditol
"""
from rdkit import Chem

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol having the general formula HOCH2[CH(OH)]nCH2OH,
    formally derivable from an aldose by reduction of the carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is acyclic
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic"

    # Get all carbon atoms
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    oxygen_counts = 0

    for atom in carbons:
        # Check if carbon is sp3 hybridized (single bonds only)
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            return False, "Carbon atom not sp3 hybridized"

        # Check for attached hydroxyl group
        has_OH = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    if neighbor.GetTotalNumHs() == 1:
                        has_OH = True
                        oxygen_counts += 1
                        break
        if not has_OH:
            return False, "Carbon atom without hydroxyl group found"

    # Check if terminal carbons are CH2OH
    terminal_carbons = [atom for atom in carbons if atom.GetDegree() == 2]
    if len(terminal_carbons) != 2:
        return False, "Incorrect number of terminal carbons"

    for atom in terminal_carbons:
        # Should be connected to two atoms: one carbon and one oxygen (OH group)
        neighbors = atom.GetNeighbors()
        if len(neighbors) != 2:
            return False, "Terminal carbon has incorrect number of neighbors"
        has_CH2OH = False
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and neighbor.GetTotalNumHs() == 1:
                    has_CH2OH = True
                    break
        if not has_CH2OH:
            return False, "Terminal carbon is not CH2OH group"

    return True, "Molecule is an alditol"