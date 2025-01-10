"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: epoxy fatty acid
"""
from rdkit import Chem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    An epoxy fatty acid is a heterocyclic fatty acid containing an epoxide ring as part of its structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, 'Invalid SMILES string'

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1,-]')
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_matches:
        return False, "No carboxylic acid group found"
    carboxyl_carbon_idx = carboxylic_matches[0][0]

    # Check for epoxide ring (three-membered cyclic ether)
    epoxide_pattern = Chem.MolFromSmarts('C1OC1')
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    if not epoxide_matches:
        return False, "No epoxide ring found"

    # Ensure that the epoxide ring is part of an aliphatic chain (not aromatic)
    is_aliphatic_epoxide = False
    for match in epoxide_matches:
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in match]
        if all(not atom.IsInRingSize(6) and not atom.GetIsAromatic() for atom in atoms_in_ring):
            is_aliphatic_epoxide = True
            epoxide_atom_indices = match
            break
    if not is_aliphatic_epoxide:
        return False, "Epoxide ring is not part of an aliphatic chain"

    # Check that the epoxide ring is connected to the hydrocarbon chain with carboxylic acid group
    # Verify that there is a path consisting mostly of carbon atoms
    is_connected = False
    for epoxide_atom_idx in epoxide_atom_indices:
        try:
            path = Chem.rdmolops.GetShortestPath(mol, carboxyl_carbon_idx, epoxide_atom_idx)
            # Check that the path is composed of carbon and oxygen atoms (excluding other heteroatoms)
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() in [6,8] for idx in path):
                is_connected = True
                break
        except:
            continue
    if not is_connected:
        return False, "Epoxide ring is not connected to the hydrocarbon chain with carboxylic acid group"

    # Count the number of carbon atoms in the chain between the epoxide and carboxylic acid
    carbon_chain_length = sum(1 for idx in path if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    if carbon_chain_length < 8:
        return False, f"Carbon chain length is {carbon_chain_length}, which is less than 8"

    # Check for aromatic rings in the molecule
    if mol.GetRingInfo().NumAromaticRings() > 0:
        return False, "Molecule contains aromatic rings"

    return True, "Molecule contains a carboxylic acid group and an epoxide ring connected by a hydrocarbon chain"

__metadata__ = {
    'chemical_class': {
        'name': 'epoxy fatty acid',
        'definition': 'A heterocyclic fatty acid containing an epoxide ring as part of its structure.'
    }
}