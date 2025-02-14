"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is an organic compound where a single sp3-hybridized carbon atom is attached to both
    an amino group and a hydroxy group, and not part of a carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise.
        str: Reason for classification.
    """

    # Parse the SMILES string to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Iterate over all carbon atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # Skip if not carbon

        if atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue  # Skip if not sp3-hybridized

        # Exclude carbons that are part of a carbonyl group (C=O)
        is_carbonyl = False
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:  # Oxygen
                    is_carbonyl = True
                    break
        if is_carbonyl:
            continue

        # Initialize counts
        has_hydroxyl = False
        has_amino = False
        other_heteroatom = False

        # Analyze neighbors
        for neighbor in atom.GetNeighbors():
            nbr_atomic_num = neighbor.GetAtomicNum()
            if nbr_atomic_num == 8:
                # Possible hydroxyl group
                if neighbor.GetTotalDegree() == 1 and neighbor.GetTotalNumHs() == 1:
                    has_hydroxyl = True
                else:
                    other_heteroatom = True  # Oxygen that's not hydroxyl
            elif nbr_atomic_num == 7:
                # Possible amino group
                if neighbor.GetFormalCharge() == 0:
                    has_amino = True
                else:
                    other_heteroatom = True  # Charged nitrogen
            elif nbr_atomic_num not in [1, 6]:  # Exclude hydrogen and carbon
                other_heteroatom = True

        # Check if the atom meets the criteria for a hemiaminal carbon
        if has_hydroxyl and has_amino and not other_heteroatom:
            return True, "Contains a carbon atom attached to both hydroxyl and amino groups (hemiaminal)."

    return False, "Does not contain a hemiaminal functional group."