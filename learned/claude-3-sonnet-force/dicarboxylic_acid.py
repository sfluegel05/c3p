"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is any carboxylic acid containing two carboxy groups
    connected by a linear or cyclic carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a dicarboxylic acid, False otherwise
        str: Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for two carboxyl groups
    carboxyl_pattern = Chem.MolFromSmarts("[C](=O)(O)")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 2:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, expected 2"

    # Check if the carboxyl groups are connected by a linear or cyclic carbon chain
    carboxyl_atoms = [mol.GetAtomWithIdx(match[0]) for match in carboxyl_matches]
    if not are_carboxyls_connected(mol, carboxyl_atoms):
        return False, "Carboxyl groups are not connected by a carbon chain"

    # Additional checks or filters (optional)
    # ...

    return True, "Carboxyl groups are connected by a carbon chain"

def are_carboxyls_connected(mol, carboxyl_atoms):
    """
    Checks if the given carboxyl atoms are connected by a linear or cyclic carbon chain.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit Mol object
        carboxyl_atoms (list): List of RDKit Atom objects representing carboxyl groups

    Returns:
        bool: True if the carboxyl groups are connected by a carbon chain, False otherwise
    """
    # Get the neighbors of each carboxyl atom
    carboxyl_neighbors = [set(atom.GetNeighbors()) for atom in carboxyl_atoms]

    # Find the common neighbors between the carboxyl atoms
    common_neighbors = carboxyl_neighbors[0].intersection(*carboxyl_neighbors[1:])

    # If there are common neighbors, the carboxyl groups are connected by a carbon chain
    if common_neighbors:
        # Check if the common neighbors form a linear or cyclic carbon chain
        chain_atoms = list(common_neighbors)
        for atom in chain_atoms:
            if not atom.GetAtomicNum() == 6:  # Carbon
                return False

        # Check for linearity or cyclicity
        if len(chain_atoms) > 2:
            # Check for linearity
            for i in range(len(chain_atoms) - 1):
                if len(set(chain_atoms[i].GetNeighbors()).intersection(set(chain_atoms[i + 1].GetNeighbors()))) != 2:
                    break
            else:
                return True
            # Check for cyclicity
            ring_info = mol.GetRingInfo()
            for ring in ring_info.AtomRings():
                if set(ring).issuperset(set(chain_atoms)):
                    return True

    return False