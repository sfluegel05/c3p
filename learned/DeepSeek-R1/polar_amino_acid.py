"""
Classifies: CHEBI:26167 polar amino acid
"""
"""
Classifies: Polar amino acids (CHEBI: ???)
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    Polar amino acids have side chains capable of forming hydrogen bonds.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Amino acid pattern: amino group (NH2/NH3+) connected to alpha carbon with carboxyl group
    amino_acid_pattern = MolFromSmarts("[NX3;H2,H3+]-[CX4H]([CX3](=[OX1])[OX2H1,O-])-[!H0]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Not an amino acid"
    
    # Find alpha carbon (second atom in the match)
    matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not matches:
        return False, "No alpha carbon found"
    alpha_carbon_idx = matches[0][1]
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
    
    # Identify side chain starting atom (non-amino, non-carboxyl substituent)
    side_chain_start = None
    for neighbor in alpha_carbon.GetNeighbors():
        # Check if neighbor is part of amino group
        if neighbor.GetAtomicNum() == 7 and neighbor.GetTotalNumHs() >= 2:
            continue
        # Check if neighbor is part of carboxyl group
        if neighbor.GetAtomicNum() == 6:
            carboxyl_oxygen = any(a.GetAtomicNum() == 8 and a.GetTotalNumHs() >=1 
                                 for a in neighbor.GetNeighbors())
            if carboxyl_oxygen:
                continue
        side_chain_start = neighbor
        break
    
    if not side_chain_start:
        return False, "No side chain found"
    
    # Collect all side chain atoms (excluding alpha carbon)
    side_chain_atoms = set()
    stack = [side_chain_start]
    while stack:
        atom = stack.pop()
        if atom.GetIdx() == alpha_carbon_idx or atom.GetIdx() in side_chain_atoms:
            continue
        side_chain_atoms.add(atom.GetIdx())
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() not in side_chain_atoms and nbr.GetIdx() != alpha_carbon_idx:
                stack.append(nbr)
    
    # Check for polar groups in side chain atoms
    polar_patterns = {
        "hydroxyl": MolFromSmarts("[OX2H]"),
        "thiol": MolFromSmarts("[SX2H]"),
        "amide": MolFromSmarts("[CX3](=O)[NX3H2]"),
        "guanidine": MolFromSmarts("[NH]C(=N)N"),
        "imidazole": MolFromSmarts("[nH]1cccn1"),
        "carboxyl": MolFromSmarts("[CX3](=O)[OX2H1,O-]"),
        "amino": MolFromSmarts("[NX3H2]"),
        "ether_oxygen": MolFromSmarts("[OX2H0][#6]"),
        "aromatic_nitrogen": MolFromSmarts("[nH]"),
    }
    
    for group_name, pattern in polar_patterns.items():
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            if any(idx in side_chain_atoms for idx in match):
                return True, f"Side chain contains {group_name} group"
    
    # Special case: aromatic rings with polar substituents (e.g., tyrosine)
    for atom_idx in side_chain_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetIsAromatic():
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() == 1:
                    return True, "Aromatic ring with hydroxyl group"
    
    return False, "No polar groups in side chain"