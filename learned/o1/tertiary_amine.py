"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: CHEBI:32876 tertiary amine
"""

from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is a compound where a nitrogen atom is sp3-hybridized,
    bonded to three hydrocarbyl groups (carbon atoms), not part of an amide, imine,
    nitrile, nitro group, or aromatic system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # Check hybridization (should be sp3)
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            # Check if nitrogen is aromatic
            if atom.GetIsAromatic():
                continue
            # Check if nitrogen has any hydrogens attached
            if atom.GetTotalNumHs(includeNeighbors=True) != 0:
                continue
            # Check if nitrogen is part of unwanted functional groups
            # Exclude amides, imines, nitriles, nitro groups
            is_amide = False
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nbr.GetAtomicNum() == 6:
                    # Check for C=O (amide) or C=N (imine)
                    for obond in nbr.GetBonds():
                        if obond.GetBondType() == Chem.rdchem.BondType.DOUBLE and obond.GetOtherAtom(nbr).GetAtomicNum() in [7,8]:
                            is_amide = True
                            break
                elif bond.GetBondType() == Chem.rdchem.BondType.TRIPLE and nbr.GetAtomicNum() == 6:
                    # Nitrile group
                    is_amide = True
                    break
                elif nbr.GetAtomicNum() == 8 and nbr.GetFormalCharge() == -1:
                    # Nitro group (N(=O)[O-])
                    is_amide = True
                    break
                if is_amide:
                    break
            if is_amide:
                continue
            # Check that nitrogen is bonded to exactly three carbon atoms
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 3:
                continue
            all_carbons = all(nbr.GetAtomicNum() == 6 for nbr in neighbors)
            if not all_carbons:
                continue
            # Passed all checks, it's a tertiary amine
            return True, "Contains a tertiary amine group"

    return False, "No tertiary amine group found"