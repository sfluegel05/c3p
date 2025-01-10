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

    # Check molecular weight - amino acids should be within specific range
    mol_weight = Descriptors.ExactMolWt(mol)
    if mol_weight < 75 or mol_weight > 205:  # Adjusted range for standard amino acids
        return False, "Molecular weight outside typical amino acid range"

    # More specific amino acid backbone pattern including zwitterionic forms
    backbone_patterns = [
        "[NX3,NX4+][CX4][CX3](=[OX1])[OX2H,OX1-]",  # Standard form
        "[NH3+][CX4][CX3](=[OX1])[O-]"  # Zwitterionic form
    ]
    
    backbone_found = False
    for pattern in backbone_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            backbone_found = True
            break
            
    if not backbone_found:
        return False, "No amino acid backbone found"

    # Check for peptide bonds to exclude peptides
    peptide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[CX4][NX3]")
    if mol.HasSubstructMatch(peptide_pattern):
        return False, "Molecule appears to be a peptide"

    # Get backbone atoms
    backbone_match = mol.GetSubstructMatch(Chem.MolFromSmarts("[NX3,NX4+][CX4][CX3](=[OX1])[OX2H,OX1-]"))
    if not backbone_match:
        return False, "Cannot identify backbone atoms"
    
    backbone_atoms = set(backbone_match)

    # Patterns for polar side chain groups
    polar_patterns = {
        "hydroxyl": ("[OX2H]", lambda m: not any(idx in backbone_atoms for idx in m)),  # Ser, Thr
        "amide": ("[NX3H2][CX3](=[OX1])", lambda m: True),  # Asn, Gln
        "carboxyl": ("[CX3](=[OX1])[OX2H,OX1-]", lambda m: not all(idx in backbone_atoms for idx in m)),  # Asp, Glu
        "basic_N": ("[NX3;H2,H1;!$(NC=O)]", lambda m: not any(idx in backbone_atoms for idx in m)),  # Lys
        "guanidino": ("[NX3][CX3](=[NX2])[NX3]", lambda m: True),  # Arg
        "imidazole": ("c1c[nH]cn1", lambda m: True),  # His
        "thiol": ("[SX2H]", lambda m: True),  # Cys
        "phenol": ("c1cc(O)ccc1", lambda m: True)  # Tyr
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

    # Count total atoms to ensure we're not dealing with modified amino acids
    atom_count = mol.GetNumAtoms()
    if atom_count > 25:  # Most amino acids have fewer atoms
        return False, "Too many atoms for a standard amino acid"

    # Construct reason string
    polar_groups_str = ", ".join(found_polar_groups)
    return True, f"Polar amino acid with {polar_groups_str} group(s) in side chain"