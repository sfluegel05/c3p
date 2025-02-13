"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem

def is_polar_amino_acid(smiles: str):
    """
    Determines if a molecule is a polar amino acid based on its SMILES string.
    A polar amino acid contains side chains capable of forming hydrogen bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polar amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized amino acid backbone pattern including isotopic versions
    aa_backbone = Chem.MolFromSmarts("N[C@@H](C)C(=O)O")  # Using chirality to ensure alpha-carbon is included
    if not mol.HasSubstructMatch(aa_backbone):
        return False, "No amino acid backbone found"
    
    # Potential polar side chain patterns including exclusions of backbone components
    polar_side_chain_patterns = [
        "[CX4][CX3](O)[OH2]",  # Alcohol/hydroxyl, e.g., serine
        "[NX3][CX3]=[OX1]",    # Amide, e.g., asparagine
        "[S][CX4]",            # Thiol, e.g., cysteine
        "[nH]",                # Imidazole, e.g., histidine
        "[OX1][CX3]=O",        # Carboxylate, e.g., aspartic acid, in context of side chain
        "[CX4][CX4][NX3]",     # Amine group in side chain, e.g., lysine
    ]

    # Check out side chain functionality
    mol_fragments = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)

    for frag in mol_fragments:
        if frag.HasSubstructMatch(aa_backbone):
            continue

        for pattern in polar_side_chain_patterns:
            polar_group = Chem.MolFromSmarts(pattern)
            if frag.HasSubstructMatch(polar_group):
                return True, f"Polar group '{pattern}' found in side chain"

    return False, "No polar groups found in side chain"