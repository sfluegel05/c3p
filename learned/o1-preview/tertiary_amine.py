"""
Classifies: CHEBI:32876 tertiary amine
"""
from rdkit import Chem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule is a tertiary amine based on its SMILES string.
    A tertiary amine is a nitrogen atom bonded to three carbon atoms via single bonds,
    with no hydrogen atoms attached, not in an aromatic ring, and zero formal charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all atoms in the molecule
    for atom in mol.GetAtoms():
        # Check if atom is nitrogen
        if atom.GetAtomicNum() == 7:
            # Check that nitrogen is not in an aromatic ring
            if not atom.GetIsAromatic():
                # Check that nitrogen has zero formal charge
                if atom.GetFormalCharge() == 0:
                    # Check that nitrogen has degree 3 (connected to 3 atoms)
                    if atom.GetDegree() == 3:
                        # Check that nitrogen has zero implicit and explicit hydrogens
                        if atom.GetTotalNumHs() == 0:
                            # Check that nitrogen is connected to three carbon atoms via single bonds
                            is_tertiary_amine = True
                            for neighbor in atom.GetNeighbors():
                                # Check that neighbor is carbon
                                if neighbor.GetAtomicNum() != 6:
                                    is_tertiary_amine = False
                                    break
                                # Check that bond is single
                                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                                if bond.GetBondType() != Chem.BondType.SINGLE:
                                    is_tertiary_amine = False
                                    break
                            if is_tertiary_amine:
                                return True, "Contains a tertiary amine nitrogen"
    return False, "No tertiary amine nitrogen found"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32877',
                          'name': 'tertiary amine',
                          'definition': 'A compound formally derived from ammonia by replacing three hydrogen atoms by hydrocarbyl groups.',
                          'parents': ['CHEBI:32874']}}