"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic amino acid structure (NH2-CH-COOH)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4][CX3](=[OX1])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid backbone found"

    # Initialize patterns for polar groups
    patterns = {
        "hydroxyl": "[OX2H]",  # -OH (serine, threonine, tyrosine)
        "amide": "[NX3H2][CX3](=[OX1])",  # -CONH2 (asparagine, glutamine)
        "carboxyl": "[CX3](=[OX1])[OX2H,OX1-]",  # -COOH (aspartic acid, glutamic acid)
        "basic_N": "[NX3;H2,H1;!$(NC=O)]",  # Basic N (lysine, arginine)
        "guanidino": "[NX3][CX3](=[NX2])[NX3]",  # Guanidino group (arginine)
        "imidazole": "c1c[nH]cn1",  # Imidazole ring (histidine)
        "thiol": "[SX2H]"  # -SH (cysteine)
    }

    # Count polar groups (excluding those in the amino acid backbone)
    found_polar_groups = []
    backbone_atoms = set()
    
    # Get backbone atoms to exclude them from side chain search
    for match in mol.GetSubstructMatches(amino_acid_pattern):
        backbone_atoms.update(match)
    
    # Check each polar pattern
    for group_name, pattern in patterns.items():
        pattern_mol = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(pattern_mol)
        
        # Only count matches that aren't part of the backbone
        side_chain_matches = [match for match in matches if not all(atom_idx in backbone_atoms for atom_idx in match)]
        
        if side_chain_matches:
            found_polar_groups.append(group_name)

    if not found_polar_groups:
        return False, "No polar side chain groups found"

    # Construct reason string
    polar_groups_str = ", ".join(found_polar_groups)
    return True, f"Polar amino acid with {polar_groups_str} group(s) in side chain"