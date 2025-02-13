"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
"""
Classifies: CHEBI:16411 branched-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens
    AllChem.EmbedMolecule(mol)
    AllChem.AddHs(mol)

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Find the longest carbon chain
    longest_chain = find_longest_carbon_chain(mol)
    if longest_chain is None or len(longest_chain) < 4:
        return False, "No carbon chain of sufficient length found"

    # Check for branches
    branches = find_branches(mol, longest_chain)
    if len(branches) < 1:
        return False, "No branches found on the carbon chain"

    # Check for cyclopropyl rings
    has_cyclopropyl_rings = check_cyclopropyl_rings(mol)

    # Check molecular weight and composition
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 150 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for fatty acids"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count < 6 or o_count < 2:
        return False, "Composition does not match a fatty acid"

    reason = "Contains a carbon chain with at least one branch"
    if has_cyclopropyl_rings:
        reason += ", and cyclopropyl ring(s)"

    return True, reason

def find_longest_carbon_chain(mol):
    """
    Finds the longest continuous carbon chain in a molecule.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object

    Returns:
        list: List of atom indices representing the longest carbon chain, or None if not found
    """
    longest_chain = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            chain = find_chain(mol, atom.GetIdx())
            if len(chain) > len(longest_chain):
                longest_chain = chain
    return longest_chain

def find_chain(mol, start_idx):
    """
    Finds a continuous chain of carbon atoms starting from a given atom index.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object
        start_idx (int): Index of the starting atom

    Returns:
        list: List of atom indices representing the chain
    """
    chain = [start_idx]
    visited = set(chain)
    to_visit = list(atom.GetNeighborIndices() for atom in mol.GetAtoms())

    while to_visit:
        neighbors = to_visit.pop(0)
        for neighbor in neighbors:
            if neighbor not in visited and mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6:
                visited.add(neighbor)
                chain.append(neighbor)
                to_visit.append(mol.GetAtomWithIdx(neighbor).GetNeighborIndices())

    return chain

def find_branches(mol, longest_chain):
    """
    Finds the branches on the longest carbon chain.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object
        longest_chain (list): List of atom indices representing the longest carbon chain

    Returns:
        list: List of atom indices representing the branches
    """
    branches = []
    for idx in longest_chain:
        atom = mol.GetAtomWithIdx(idx)
        neighbors = [mol.GetAtomWithIdx(n_idx) for n_idx in atom.GetNeighborIndices()]
        carbon_neighbors = [n for n in neighbors if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) > 2:
            for n in carbon_neighbors:
                if n.GetIdx() not in longest_chain:
                    branches.append(n.GetIdx())
    return list(set(branches))

def check_cyclopropyl_rings(mol):
    """
    Checks if the molecule contains cyclopropyl rings.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object

    Returns:
        bool: True if the molecule contains cyclopropyl rings, False otherwise
    """
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        carbon_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 6]
        if len(carbon_atoms) == 3:
            return True
    return False