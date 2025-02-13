"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: prostaglandin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    A prostaglandin is a naturally occurring compound derived from prostanoic acid,
    characterized by a C20 backbone with a cyclopentane ring and two side chains,
    one of which ends with a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check total number of carbons (should be approximately 20)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 17 or num_carbons > 24:
        return False, f"Number of carbons is {num_carbons}, which is not in the expected range for prostaglandins"

    # Get ring information
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return False, f"Molecule has {ring_info.NumRings()} rings, expected 1 cyclopentane ring"

    # Check if the ring is a 5-membered ring
    atom_rings = ring_info.AtomRings()
    ring_sizes = [len(ring) for ring in atom_rings]
    if ring_sizes[0] != 5:
        return False, f"Ring size is {ring_sizes[0]}, expected 5"

    # Check for carboxylic acid group
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[OH]')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Check for hydroxyl group(s)
    hydroxyl = Chem.MolFromSmarts('[OX2H]')
    if not mol.HasSubstructMatch(hydroxyl):
        return False, "No hydroxyl groups found"

    # Check that the cyclopentane ring has at least two side chains attached
    cp_ring_atoms = atom_rings[0]
    side_chain_attachments = []
    for idx in cp_ring_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in cp_ring_atoms:
                side_chain_attachments.append(idx)
                break  # Only need to know if there is at least one neighbor outside the ring

    if len(side_chain_attachments) < 2:
        return False, "Less than two side chains attached to cyclopentane ring"

    # If all checks pass, it is likely a prostaglandin
    return True, "Molecule matches prostaglandin structural features"

__metadata__ = {
    'chemical_class': {
        'name': 'prostaglandin',
        'definition': 'Naturally occurring compounds derived from the parent C20 acid, prostanoic acid.',
        'parents': []
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
    },
    'message': None,
    'attempt': 2,
    'success': True,
    'error': '',
    'stdout': None
}