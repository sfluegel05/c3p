"""
Classifies: CHEBI:26167 polar amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

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

    # Check molecular weight - amino acids should be relatively small
    mol_weight = Descriptors.ExactMolWt(mol)
    if mol_weight > 250:  # Most amino acids are under 200 Da
        return False, "Molecule too large to be a single amino acid"

    # Check for basic amino acid structure (NH2-CH-COOH)
    amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4][CX3](=[OX1])[OX2H,OX1-]")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid backbone found"

    # Check number of amino groups and carboxyl groups
    amino_pattern = Chem.MolFromSmarts("[NX3,NX4+][CX4][CX3](=[OX1])[OX2H,OX1-]")
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX1-]")
    
    amino_count = len(mol.GetSubstructMatches(amino_pattern))
    carboxyl_count = len(mol.GetSubstructMatches(carboxyl_pattern))
    
    if amino_count > 1 or carboxyl_count > 2:  # Allow up to 2 carboxyls for asp/glu
        return False, "Multiple amino acid residues detected - likely a peptide"

    # Initialize patterns for polar groups in side chains
    patterns = {
        "hydroxyl": "[OX2H]",  # -OH (serine, threonine, tyrosine)
        "amide": "[NX3H2][CX3](=[OX1])",  # -CONH2 (asparagine, glutamine)
        "carboxyl": "[CX3](=[OX1])[OX2H,OX1-]",  # -COOH (aspartic acid, glutamic acid)
        "basic_N": "[NX3;H2,H1;!$(NC=O)]",  # Basic N (lysine, arginine)
        "guanidino": "[NX3][CX3](=[NX2])[NX3]",  # Guanidino group (arginine)
        "imidazole": "c1c[nH]cn1",  # Imidazole ring (histidine)
        "thiol": "[SX2H]",  # -SH (cysteine)
        "indole": "c1ccc2c(c1)[nH]cc2",  # Indole ring (tryptophan)
        "phenol": "c1cc(O)ccc1"  # Phenol group (tyrosine)
    }

    # Get backbone atoms to exclude them from side chain search
    backbone_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not backbone_matches:
        return False, "No amino acid backbone found"
    
    backbone_atoms = set(backbone_matches[0])
    
    # Check each polar pattern
    found_polar_groups = []
    for group_name, pattern in patterns.items():
        pattern_mol = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(pattern_mol)
        
        # Only count matches that aren't part of the backbone
        side_chain_matches = [match for match in matches if not all(atom_idx in backbone_atoms for atom_idx in match)]
        
        if side_chain_matches:
            found_polar_groups.append(group_name)

    if not found_polar_groups:
        # Special case for tryptophan - check for indole NH
        if mol.HasSubstructMatch(Chem.MolFromSmarts("c1ccc2c(c1)[nH]cc2")):
            return True, "Polar amino acid with indole NH group capable of hydrogen bonding"
        return False, "No polar side chain groups found"

    # Construct reason string
    polar_groups_str = ", ".join(found_polar_groups)
    return True, f"Polar amino acid with {polar_groups_str} group(s) in side chain"