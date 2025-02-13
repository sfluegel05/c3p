"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:35612 secondary amine
A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of nitrogen atom
    if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7) == 0:
        return False, "No nitrogen atom found"

    # Check for single nitrogen atom bonded to exactly two carbons
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    is_secondary_amine = False
    for n_atom in nitrogen_atoms:
        n_neighbors = [neighbor for neighbor in n_atom.GetNeighbors() if neighbor.GetAtomicNum() == 6]
        if len(n_neighbors) == 2:
            is_secondary_amine = True
            break

    if not is_secondary_amine:
        return False, "Nitrogen atom not bonded to exactly two carbons"

    # Check for two separate alkyl/aryl groups attached to nitrogen
    n_atom = n_neighbors[0]
    first_chain = [n_atom]
    second_chain = [n_neighbors[1]]

    while True:
        next_atoms = []
        for atom in first_chain:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor not in first_chain:
                    next_atoms.append(neighbor)
        if not next_atoms:
            break
        first_chain.extend(next_atoms)

    while True:
        next_atoms = []
        for atom in second_chain:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor not in second_chain:
                    next_atoms.append(neighbor)
        if not next_atoms:
            break
        second_chain.extend(next_atoms)

    if len(first_chain) == 1 or len(second_chain) == 1:
        return False, "At least one of the substituents is not an alkyl/aryl group"

    return True, "Contains a nitrogen atom bonded to two separate alkyl/aryl groups"