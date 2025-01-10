"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define quinone patterns (1,4-benzoquinone and 1,4-naphthoquinone)
    quinone_patterns = [
        Chem.MolFromSmarts('O=C1C=CC=CC1=O'),  # 1,4-benzoquinone
        Chem.MolFromSmarts('O=C1C=CC2=CC=CC=C2C1=O')  # 1,4-naphthoquinone
    ]
    
    # Check for presence of quinone ring
    quinone_match = False
    for pattern in quinone_patterns:
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            quinone_match = True
            break
    if not quinone_match:
        return False, "No quinone ring found"

    # Get ring atoms from the match
    ring_atoms = set()
    for idx in matches[0]:
        ring_atoms.add(idx)

    # Find side chains attached to the quinone ring
    side_chains = []
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in ring_atoms:
                # Found substituent attached to ring
                # Traverse the side chain starting from this neighbor atom
                # Exclude ring atoms
                atoms_in_side_chain = set()
                bonds_in_side_chain = set()
                stack = [neighbor_idx]
                while stack:
                    current_idx = stack.pop()
                    if current_idx not in atoms_in_side_chain:
                        atoms_in_side_chain.add(current_idx)
                        current_atom = mol.GetAtomWithIdx(current_idx)
                        for nbr in current_atom.GetNeighbors():
                            nbr_idx = nbr.GetIdx()
                            bond = mol.GetBondBetweenAtoms(current_idx, nbr_idx)
                            if nbr_idx not in atoms_in_side_chain and nbr_idx not in ring_atoms:
                                stack.append(nbr_idx)
                                bonds_in_side_chain.add(bond.GetIdx())
                # Analyze the side chain
                num_carbons = 0
                num_double_bonds = 0
                for idx in atoms_in_side_chain:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() == 6:
                        num_carbons +=1
                for bond_idx in bonds_in_side_chain:
                    bond = mol.GetBondWithIdx(bond_idx)
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        num_double_bonds +=1
                side_chains.append({
                    'num_carbons': num_carbons,
                    'num_double_bonds': num_double_bonds
                })

    # Check if any side chain meets the criteria for a polyprenyl chain
    for chain in side_chains:
        if chain['num_carbons'] >= 10 and chain['num_double_bonds'] >= 2:
            return True, "Contains quinone ring with polyprenyl side chain"

    return False, "No polyprenyl side chain attached to quinone ring"