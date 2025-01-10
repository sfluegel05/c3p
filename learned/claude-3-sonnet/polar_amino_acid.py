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

    # Check for modifications that would exclude the molecule
    exclusion_patterns = [
        ("[NX3](C(=O)[CH3])", "N-acetylated amino acid"),  # N-acetyl
        ("[CX4]([NX3])([CX3](=O)[OH])[CX4][OH]", "alpha-hydroxy amino acid"),  # alpha-hydroxy
        ("[B,Si,P]", "Contains B, Si, or P atoms"),  # Non-standard elements
    ]
    
    for pattern, reason in exclusion_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, f"Not a standard amino acid: {reason}"

    # Strict amino acid backbone pattern - ensures NH2/NH3+ is directly on alpha carbon
    backbone_patterns = [
        "[NX3H2,NX4H3+][CX4H1]([*])[CX3](=[OX1])[OX2H,OX1-]",  # Standard/zwitterionic
    ]
    
    backbone_found = False
    backbone_match = None
    for pattern in backbone_patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(pattern_mol)
        if matches:
            backbone_found = True
            backbone_match = matches[0]
            break
            
    if not backbone_found:
        return False, "No valid amino acid backbone found"

    # Get backbone atoms including the alpha carbon substituent
    backbone_atoms = set(backbone_match)
    alpha_carbon_idx = backbone_match[1]  # Index of alpha carbon from backbone match

    # Patterns for polar side chain groups with improved specificity
    polar_patterns = {
        "hydroxyl": ("[OX2H][CX4]", lambda m: not any(idx in backbone_atoms for idx in m)),  # Ser, Thr
        "amide": ("[NX3H2][CX3](=[OX1])", lambda m: not any(idx in backbone_atoms for idx in m)),  # Asn, Gln
        "carboxyl": ("[CX3](=[OX1])[OX2H,OX1-]", lambda m: not any(idx in backbone_atoms for idx in m)),  # Asp, Glu
        "basic_N": ("[NX3H2][CX4]", lambda m: not any(idx in backbone_atoms for idx in m)),  # Lys
        "guanidino": ("[NX3][CX3](=[NX2])[NX3]", lambda m: True),  # Arg
        "imidazole": ("c1c[nH]cn1", lambda m: True),  # His
        "thiol": ("[SX2H]", lambda m: not any(idx in backbone_atoms for idx in m)),  # Cys
        "phenol": ("c1cc(O)ccc1", lambda m: True),  # Tyr
        "indole": ("c1ccc2c(c1)[nH]cc2", lambda m: True)  # Trp
    }

    # Check for polar groups in side chain
    found_polar_groups = []
    for group_name, (pattern, validator) in polar_patterns.items():
        pattern_mol = Chem.MolFromSmarts(pattern)
        matches = mol.GetSubstructMatches(pattern_mol)
        
        # Validate matches using the custom validator function
        valid_matches = [m for m in matches if validator(m)]
        
        if valid_matches:
            found_polar_groups.append(group_name)

    if not found_polar_groups:
        return False, "No polar side chain groups found"

    # Additional validation for molecular properties
    if mol.GetNumAtoms() > 30:  # Stricter atom count limit
        return False, "Too many atoms for a standard amino acid"
    
    mol_weight = Descriptors.ExactMolWt(mol)
    if mol_weight < 75 or mol_weight > 250:  # Adjusted range
        return False, "Molecular weight outside typical amino acid range"

    # Construct reason string
    polar_groups_str = ", ".join(found_polar_groups)
    return True, f"Polar amino acid with {polar_groups_str} group(s) in side chain"