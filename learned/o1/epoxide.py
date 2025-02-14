"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:32987 epoxide
"""
from rdkit import Chem

def is_epoxide(smiles: str):
    """
    Determines if a molecule is an epoxide based on its SMILES string.
    An epoxide is any cyclic ether in which the oxygen atom forms part of a 3-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxide, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all rings in the molecule
    sssr = Chem.GetSymmSSSR(mol)

    # Flag to indicate if an epoxide ring is found
    epoxide_found = False

    # Iterate over each ring
    for ring in sssr:
        # Check if the ring has 3 atoms
        if len(ring) == 3:
            # Count the number of oxygen atoms in the ring
            oxygen_count = 0
            carbon_count = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:  # Oxygen
                    oxygen_count += 1
                elif atom.GetAtomicNum() == 6:  # Carbon
                    carbon_count += 1
            # Check if the ring has exactly one oxygen and two carbons
            if oxygen_count == 1 and carbon_count == 2:
                epoxide_found = True
                break

    if epoxide_found:
        return True, "Contains an epoxide ring (3-membered ring with one oxygen atom)"
    else:
        return False, "Does not contain an epoxide ring"