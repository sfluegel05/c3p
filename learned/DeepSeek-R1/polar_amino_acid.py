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
    
    # Check for amino acid structure: alpha-amino and alpha-carboxyl groups
    amino_acid_pattern = MolFromSmarts("[NH2][C]([CX3]=O)([!H])[!H]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "Not an amino acid"
    
    # Identify the side chain (R group) attached to the alpha carbon
    alpha_carbon = MolFromSmarts("[NH2][C]([CX3]=O)([!H])[R]")
    matches = mol.GetSubstructMatches(alpha_carbon)
    if not matches:
        return False, "No side chain found"
    
    # Extract the side chain atoms
    side_chain_atoms = []
    for match in matches:
        # The fourth atom in the match is the R group attachment point
        r_attachment = match[3]
        # Traverse the side chain starting from the attachment
        visited = set()
        stack = [r_attachment]
        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in visited and neighbor.GetIdx() not in [match[0], match[1], match[2]]:
                    stack.append(neighbor.GetIdx())
        side_chain_atoms = list(visited)
        break  # Assuming one alpha carbon
    
    # Create a side chain molecule
    side_chain = Chem.PathToSubmol(mol, side_chain_atoms)
    if not side_chain:
        return False, "No side chain"
    
    # Check for polar groups in the side chain
    polar_patterns = {
        "hydroxyl": MolFromSmarts("[OH]"),
        "thiol": MolFromSmarts("[SH]"),
        "amide": MolFromSmarts("[CX3](=O)[NX3H2]"),
        "guanidine": MolFromSmarts("[NH]C(=N)N"),
        "imidazole": MolFromSmarts("c1cncn1"),
        "carboxyl": MolFromSmarts("[CX3](=O)[OX2H1]"),
        "amino": MolFromSmarts("[NH2]")
    }
    
    for group, pattern in polar_patterns.items():
        if side_chain.HasSubstructMatch(pattern):
            return True, f"Side chain contains {group} group"
    
    # Check for aromatic rings with polar substituents (e.g., tyrosine)
    aromatic_rings = Chem.GetSymmSSSR(side_chain)
    for ring in aromatic_rings:
        for atom_idx in ring:
            atom = side_chain.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 8 or atom.GetAtomicNum() == 7:
                return True, "Aromatic ring with polar substituent"
    
    return False, "No polar groups in side chain"