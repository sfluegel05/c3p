"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone based on its SMILES string.
    A secondary alpha-hydroxy ketone has a hydroxy group attached to the alpha carbon adjacent to a ketone group,
    where the alpha carbon is secondary (attached to two carbons and one hydrogen).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all oxygen atoms to find hydroxyl groups
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen atom
            # Check if oxygen is part of a hydroxyl group (single bond to carbon)
            if atom.GetDegree() == 1:
                neighbor = atom.GetNeighbors()[0]
                if neighbor.GetAtomicNum() == 6:
                    alpha_carbon = neighbor
                    # Check bond type between oxygen and carbon
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), alpha_carbon.GetIdx())
                    if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                        continue

                    # Check if alpha carbon is secondary (attached to two carbons and one hydrogen)
                    num_hydrogens = alpha_carbon.GetTotalNumHs()
                    num_carbon_neighbors = sum(1 for nbr in alpha_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6)
                    if num_hydrogens != 1 or num_carbon_neighbors != 2:
                        continue

                    # Check if alpha carbon is adjacent to a ketone group (C=O)
                    found_ketone = False
                    for nbr in alpha_carbon.GetNeighbors():
                        if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != neighbor.GetIdx():
                            for bond in nbr.GetBonds():
                                other_atom = bond.GetOtherAtom(nbr)
                                if (other_atom.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE):
                                    found_ketone = True
                                    break
                            if found_ketone:
                                break

                    if found_ketone:
                        return True, "Contains secondary alpha-hydroxy ketone substructure"

    return False, "Does not contain secondary alpha-hydroxy ketone substructure"