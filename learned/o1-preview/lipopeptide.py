"""
Classifies: CHEBI:46895 lipopeptide
"""
"""
Classifies: CHEBI:25212 lipopeptide

"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopeptide(smiles: str):
    """
    Determines if a molecule is a lipopeptide based on its SMILES string.
    A lipopeptide is a compound consisting of a peptide with attached lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipopeptide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify peptide bonds (amide bonds connecting amino acids)
    peptide_bond_pattern = Chem.MolFromSmarts("N[C](=O)[C]")  # N-C(=O)-C
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    if len(peptide_bonds) < 2:
        return False, "Not enough peptide bonds found (need at least 2)"
    
    # Function to find the longest aliphatic carbon chain
    def get_longest_aliphatic_chain_length(mol):
        max_chain_len = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6 and not atom.IsInRing():
                visited = set()
                queue = [(atom.GetIdx(), 1)]
                while queue:
                    idx, length = queue.pop(0)
                    if idx in visited:
                        continue
                    visited.add(idx)
                    if length > max_chain_len:
                        max_chain_len = length
                    atom_idx = mol.GetAtomWithIdx(idx)
                    for neighbor in atom_idx.GetNeighbors():
                        nidx = neighbor.GetIdx()
                        n_atom = mol.GetAtomWithIdx(nidx)
                        if n_atom.GetAtomicNum() == 6 and not n_atom.IsInRing():
                            queue.append((nidx, length+1))
        return max_chain_len

    max_chain_len = get_longest_aliphatic_chain_length(mol)
    if max_chain_len < 8:
        return False, f"No long aliphatic chain found (longest chain length {max_chain_len})"

    # Check if lipid chain is connected to the peptide
    # Get atom indices for peptide bonds
    peptide_atoms = set()
    for match in peptide_bonds:
        peptide_atoms.update(match)
    
    # Get all aliphatic carbons in the longest chain
    lipid_atoms = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and not atom.IsInRing():
            lipid_atoms.add(atom.GetIdx())
    
    # For each lipid atom, check if connected to peptide
    paths_exist = False
    for lipid_atom in lipid_atoms:
        for peptide_atom in peptide_atoms:
            try:
                path = Chem.rdmolops.GetShortestPath(mol, lipid_atom, peptide_atom)
                if path:
                    paths_exist = True
                    break
            except:
                continue
        if paths_exist:
            break

    if not paths_exist:
        return False, "Lipid chain not connected to peptide moiety"

    return True, "Contains peptide bonds and attached lipid chain"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:25212',
        'name': 'lipopeptide',
        'definition': 'A compound consisting of a peptide with attached lipid.',
        'parents': []
    }
}