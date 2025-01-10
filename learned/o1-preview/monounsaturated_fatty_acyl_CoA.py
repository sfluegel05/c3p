"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
"""
Classifies: monounsaturated fatty acyl-CoA
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Check for adenine ring (as part of CoA)
    adenine_pattern = Chem.MolFromSmarts('n1cnc2c1ncnc2')
    if not mol.HasSubstructMatch(adenine_pattern):
        return False, 'No adenine moiety found (not a CoA derivative)'

    # Find thioester linkage [C(=O)S]
    thioester_pattern = Chem.MolFromSmarts('C(=O)S')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, 'No thioester linkage found'

    # Assume the fatty acyl chain is attached via the thioester linkage to CoA
    for match in thioester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon index
        sulfur_idx = match[2]      # Sulfur atom index

        # From the carbonyl carbon, collect the acyl chain atoms, excluding the sulfur atom
        acyl_chain_atom_indices = collect_acyl_chain(mol, carbonyl_c_idx, sulfur_idx)

        if not acyl_chain_atom_indices:
            continue  # Try the next thioester linkage

        # Create the acyl chain sub-molecule
        acyl_chain_mol = Chem.PathToSubmol(mol, acyl_chain_atom_indices)
        
        # Initialize ring information to prevent RingInfo error
        Chem.GetSSSR(acyl_chain_mol)

        # Analyze the acyl chain
        if acyl_chain_mol.GetRingInfo().NumRings() > 0:
            # Acyl chain contains rings
            continue

        # Count number of carbon atoms (excluding the carbonyl carbon)
        carbon_atoms = [atom for atom in acyl_chain_mol.GetAtoms() if atom.GetAtomicNum() == 6]
        num_carbons = len(carbon_atoms) - 1  # Exclude carbonyl carbon
        if num_carbons < 4:
            # Acyl chain too short
            continue

        # Count carbon-carbon double bonds in the acyl chain
        num_double_bonds = 0
        for bond in acyl_chain_mol.GetBonds():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    num_double_bonds += 1

        if num_double_bonds == 1:
            return True, 'Molecule is a monounsaturated fatty acyl-CoA (contains exactly one carbon-carbon double bond in the acyl chain)'
        else:
            continue  # Try the next thioester linkage

    return False, 'Could not find an acyl chain with exactly one carbon-carbon double bond'

def collect_acyl_chain(mol, start_idx, exclude_idx):
    '''
    Collects the atom indices of the acyl chain starting from the carbonyl carbon,
    excluding the sulfur atom and CoA moiety.

    Args:
        mol (Chem.Mol): RDKit molecule
        start_idx (int): Starting atom index (carbonyl carbon)
        exclude_idx (int): Atom index to exclude (sulfur atom)
    
    Returns:
        list: Atom indices of the acyl chain
    '''
    visited = set()
    to_visit = [start_idx]
    
    while to_visit:
        current_idx = to_visit.pop()
        if current_idx in visited:
            continue
        visited.add(current_idx)

        atom = mol.GetAtomWithIdx(current_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx == exclude_idx:
                continue
            if neighbor_idx not in visited:
                to_visit.append(neighbor_idx)

    return list(visited)

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'monounsaturated fatty acyl-CoA',
        'definition': 'Any unsaturated fatty acyl-CoA in which the fatty acyl chain contains one carbon-carbon double bond.'
    }
}