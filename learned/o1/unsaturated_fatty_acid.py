"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: unsaturated fatty acid
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid is a molecule with a long unbranched aliphatic chain attached to a carboxylic acid group,
    containing at least one carbon-carbon double or triple bond in the chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)[O;H1]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"
    
    # Consider each carboxylic acid group found
    for match in carboxy_matches:
        carboxyl_carbon_idx = match[0]
        # Find the carbon chain attached to the carboxyl carbon
        # Get the neighboring atom(s) connected to the carboxyl carbon, excluding the oxygen atoms
        carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_carbon_idx)
        neighbors = [atom for atom in carboxyl_carbon.GetNeighbors() if atom.GetAtomicNum() == 6]
        if not neighbors:
            continue  # No carbon chain attached

        chain_atom = neighbors[0]
        # Traverse the chain from this atom
        visited = set()
        chain_atoms = []

        def traverse_chain(atom, previous_atom_idx=None):
            if atom.GetIdx() in visited:
                return
            visited.add(atom.GetIdx())
            chain_atoms.append(atom)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != previous_atom_idx:
                    if not neighbor.IsInRing():
                        traverse_chain(neighbor, atom.GetIdx())
        traverse_chain(chain_atom, carboxyl_carbon_idx)
        
        # Count the number of carbons in the chain
        chain_length = len(chain_atoms) + 1  # Include the carboxyl carbon
        if chain_length < 4:
            continue  # Chain too short to be a fatty acid

        # Check for branching (atoms connected to chain atoms not in chain)
        branching = False
        for atom in chain_atoms:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor not in chain_atoms:
                    branching = True
                    break
            if branching:
                break
        if branching:
            continue  # Chain is branched, not a fatty acid

        # Check for rings in the chain
        if any(atom.IsInRing() for atom in chain_atoms):
            continue  # Chain contains ring, not a fatty acid

        # Check for C=C or C#C bonds in the chain
        has_double_triple_bond = False
        for bond in mol.GetBonds():
            if bond.GetBondType() in [Chem.rdchem.BondType.DOUBLE, Chem.rdchem.BondType.TRIPLE]:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                if begin_atom in chain_atoms and end_atom in chain_atoms:
                    has_double_triple_bond = True
                    break
        if not has_double_triple_bond:
            continue  # No carbon-carbon double or triple bond in the chain

        # If all checks passed
        return True, "Molecule is an unsaturated fatty acid"

    # If no carboxylic acid group with appropriate chain found
    return False, "Does not meet criteria for unsaturated fatty acid"