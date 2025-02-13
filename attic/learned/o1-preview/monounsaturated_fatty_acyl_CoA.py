"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: monounsaturated fatty acyl-CoA
"""

from rdkit import Chem

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    '''
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA is a fatty acyl-CoA in which the fatty acyl chain contains exactly one carbon-carbon double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    '''

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, 'Invalid SMILES string'

    # Check for thioester linkage [C(=O)-S]
    thioester_pattern = Chem.MolFromSmarts('[CX3](=O)[SX2]')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, 'No thioester linkage found'

    if len(thioester_matches) != 1:
        return False, f'Expected exactly one thioester linkage, found {len(thioester_matches)}'

    # Check for adenine ring (as part of CoA)
    adenine_pattern = Chem.MolFromSmarts('n1cnc2c1ncnc2')
    adenine_matches = mol.GetSubstructMatches(adenine_pattern)
    if not adenine_matches:
        return False, 'No adenine moiety found (not a CoA derivative)'

    # Proceed with the thioester linkage
    match = thioester_matches[0]
    carbonyl_c_idx = match[0]  # Carbonyl carbon index
    sulfur_idx = match[1]      # Sulfur index

    carbonyl_c = mol.GetAtomWithIdx(carbonyl_c_idx)

    # Find the acyl chain starting atom (connected to carbonyl carbon and not sulfur or oxygen)
    acyl_chain_start_atom = None
    for neighbor in carbonyl_c.GetNeighbors():
        if neighbor.GetAtomicNum() == 8:  # Oxygen
            continue
        if neighbor.GetAtomicNum() == 16:  # Sulfur
            continue
        # Assume this is the acyl chain starting atom
        acyl_chain_start_atom = neighbor
        break
    if acyl_chain_start_atom is None:
        return False, 'Could not find acyl chain starting atom'

    # Collect acyl chain atoms
    acyl_chain_atoms = set()
    visited_atoms = set()
    def collect_acyl_chain_atoms(atom, visited_atoms):
        visited_atoms.add(atom.GetIdx())
        acyl_chain_atoms.add(atom.GetIdx())
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() == carbonyl_c_idx:
                continue
            if neighbor.GetIdx() in visited_atoms:
                continue
            collect_acyl_chain_atoms(neighbor, visited_atoms)

    collect_acyl_chain_atoms(acyl_chain_start_atom, visited_atoms)

    # Count number of carbon-carbon double bonds in the acyl chain
    num_double_bonds = 0
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        if begin_idx in acyl_chain_atoms and end_idx in acyl_chain_atoms:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                if bond.GetBondTypeAsDouble() == 2.0:
                    num_double_bonds += 1

    if num_double_bonds == 1:
        return True, 'Molecule is a monounsaturated fatty acyl-CoA (contains exactly one carbon-carbon double bond in the acyl chain)'
    else:
        return False, f'Molecule is not a monounsaturated fatty acyl-CoA (acyl chain contains {num_double_bonds} carbon-carbon double bonds)'

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'monounsaturated fatty acyl-CoA',
        'definition': 'Any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon-carbon double bond.'
    }
}