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
    
    # Improved amino acid pattern: alpha carbon connected to amino and carboxyl groups
    # Allows for NH2, NH3+, COOH, or COO-
    amino_acid_pattern = MolFromSmarts("[NX3;H2,H3+][CX4H]([CX3]=[OX1])[!$([#1])]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Not an amino acid"
    
    # Find the alpha carbon and its substituents
    matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not matches:
        return False, "No alpha carbon found"
    
    # Get side chain atoms (the non-amino/non-carboxyl substituent on alpha carbon)
    alpha_carbon_idx = matches[0][1]
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
    
    # Find the side chain atom (non-amino, non-carboxyl neighbor of alpha carbon)
    side_chain_atom = None
    for neighbor in alpha_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 7 and neighbor.GetTotalNumHs() >= 2:
            continue  # amino group
        if neighbor.GetAtomicNum() == 6 and any(a.GetAtomicNum() == 8 for a in neighbor.GetNeighbors()):
            continue  # carboxyl group (C connected to O)
        side_chain_atom = neighbor
        break
    
    if not side_chain_atom:
        return False, "No side chain found"
    
    # Extract side chain submol starting from the side chain atom
    side_chain = Chem.ReplaceSidechains(mol, [(alpha_carbon_idx, side_chain_atom.GetIdx())])
    if not side_chain:
        return False, "Failed to extract side chain"
    
    # Check for polar groups in the side chain
    polar_patterns = {
        "hydroxyl": MolFromSmarts("[OX2H]"),  # -OH
        "thiol": MolFromSmarts("[SX2H]"),     # -SH
        "amide": MolFromSmarts("[CX3](=O)[NX3H2]"),  # -CONH2
        "guanidine": MolFromSmarts("[NH]C(=N)N"),    # arginine's group
        "imidazole": MolFromSmarts("[nH]1cccn1"),    # histidine
        "carboxyl": MolFromSmarts("[CX3](=O)[OX2H1]"),  # -COOH (as in aspartic/glutamic)
        "amino": MolFromSmarts("[NX3H2]"),           # -NH2 (lysine)
        "ether": MolFromSmarts("[OX2H0]([#6])[#6]"), # ether oxygen (e.g., serine's -CH2OH)
        "aromatic_nitrogen": MolFromSmarts("[nH]")   # aromatic N (histidine)
    }
    
    for group_name, pattern in polar_patterns.items():
        if side_chain.HasSubstructMatch(pattern):
            return True, f"Side chain contains {group_name} group"
    
    # Check for aromatic rings with polar substituents (e.g., tyrosine's -OH on benzene)
    aromatic_rings = Chem.GetSymmSSSR(side_chain)
    for ring in aromatic_rings:
        for atom_idx in ring:
            atom = side_chain.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() in [7, 8]:
                return True, "Aromatic ring with polar substituent"
    
    return False, "No polar groups in side chain"