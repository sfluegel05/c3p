"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: CHEBI:47909 2,5-diketopiperazine
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine has a piperazine-2,5-dione skeleton, which is a six-membered ring 
    with two nitrogen atoms at positions 1 and 4, and two carbonyl groups at positions 2 and 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Iterate over all 6-membered rings
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # Skip rings that are not 6-membered

        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]

        # Count nitrogen atoms in the ring
        num_nitrogens = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)

        # Count carbon atoms double-bonded to oxygen (carbonyl groups) in the ring
        num_carbonyl_carbons = 0
        for atom in ring_atoms:
            if atom.GetAtomicNum() == 6:
                # Check if this carbon atom has a double bond to oxygen
                is_carbonyl = False
                for bond in atom.GetBonds():
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        is_carbonyl = True
                        break
                if is_carbonyl:
                    num_carbonyl_carbons += 1

        # Check if ring matches 2,5-diketopiperazine skeleton
        if num_nitrogens == 2 and num_carbonyl_carbons == 2:
            return True, "Contains 2,5-diketopiperazine skeleton"

    return False, "Does not contain 2,5-diketopiperazine skeleton"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:47909',
        'name': '2,5-diketopiperazine',
        'definition': 'Any piperazinone that has a piperazine-2,5-dione skeleton.',
        'parents': ['CHEBI:24164', 'CHEBI:48373']
    },
    'config': {
        # Configuration parameters can be added here if needed
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Performance metrics can be added here if available
}