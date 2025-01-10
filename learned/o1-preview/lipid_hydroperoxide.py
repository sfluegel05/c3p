"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
"""
Classifies: lipid hydroperoxide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is any lipid carrying one or more hydroperoxy (-OOH) substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens to the molecule
    mol = Chem.AddHs(mol)

    # Function to count hydroperoxy groups, including deprotonated forms
    def count_hydroperoxy_groups(mol):
        num_hydroperoxy = 0
        for bond in mol.GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                atom1 = bond.GetBeginAtom()
                atom2 = bond.GetEndAtom()
                if atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 8:
                    # O-O single bond found
                    oxy1_h = atom1.GetTotalNumHs()
                    oxy2_h = atom2.GetTotalNumHs()
                    oxy1_charge = atom1.GetFormalCharge()
                    oxy2_charge = atom2.GetFormalCharge()
                    # Check if one oxygen has hydrogen or negative charge
                    if (oxy1_h > 0 or oxy1_charge == -1) or (oxy2_h > 0 or oxy2_charge == -1):
                        num_hydroperoxy += 1
        return num_hydroperoxy

    num_hydroperoxy = count_hydroperoxy_groups(mol)
    if num_hydroperoxy == 0:
        return False, "No hydroperoxy (-OOH) groups found"

    # Determine the length of the longest carbon chain
    def get_longest_carbon_chain_length(mol):
        # Create adjacency list of carbon atoms
        adj_list = {}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:
                idx = atom.GetIdx()
                neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
                adj_list[idx] = neighbors

        # Function to perform DFS and find longest path from a given node
        def dfs(node, visited):
            visited.add(node)
            lengths = [0]
            for neighbor in adj_list.get(node, []):
                if neighbor not in visited:
                    lengths.append(1 + dfs(neighbor, visited.copy()))
            return max(lengths)

        # Find the maximum length from each carbon atom
        max_length = 0
        for atom_idx in adj_list.keys():
            length = dfs(atom_idx, set())
            if length > max_length:
                max_length = length
        return max_length + 1  # Add 1 to include the starting atom

    longest_chain_length = get_longest_carbon_chain_length(mol)
    if longest_chain_length < 8:
        return False, f"Longest carbon chain length is {longest_chain_length}, less than 8"

    return True, f"Contains {num_hydroperoxy} hydroperoxy group(s) and a long hydrocarbon chain of length {longest_chain_length}"

__metadata__ = {
    'chemical_class': {
        'name': 'lipid hydroperoxide',
        'definition': 'Any lipid carrying one or more hydroperoxy substituents.'
    }
}