"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: essential fatty acids
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    Essential fatty acids are polyunsaturated fatty acids that must be obtained through diet.
    Common examples include omega-3 and omega-6 fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for modified compounds (deuterated, fluorinated, etc.)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8, 15, 7]:  # Allow only H,C,O,P,N
            return False, f"Contains non-standard atoms"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:  # Minimum C16 for essential fatty acids
        return False, f"Carbon chain too short ({c_count} carbons)"
    if c_count > 40:  # Maximum reasonable length
        return False, f"Carbon chain too long ({c_count} carbons)"

    # Count double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    n_double_bonds = len(double_bond_matches)
    
    if n_double_bonds < 2:  # Must be polyunsaturated
        return False, f"Not polyunsaturated (needs at least 2 double bonds)"
    if n_double_bonds > 6:
        return False, f"Too many double bonds ({n_double_bonds})"

    # Check for carboxylic acid group or valid ester in phospholipid
    carboxylic_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    phospholipid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CH2][CH]([OX2])[CH2][OX2]P(=[OX1])([OX2])")
    
    is_acid = mol.HasSubstructMatch(carboxylic_pattern)
    is_phospholipid = mol.HasSubstructMatch(phospholipid_pattern)
    
    if not (is_acid or is_phospholipid):
        return False, "Neither free fatty acid nor valid phospholipid"

    # Essential fatty acid patterns
    patterns = {
        "linoleic": Chem.MolFromSmarts("CCCCC/C=C/C/C=C/CCCCC"),  # omega-6 (9,12)
        "alpha_linolenic": Chem.MolFromSmarts("CC/C=C/C/C=C/C/C=C/CCCC"),  # omega-3 (9,12,15)
        "arachidonic": Chem.MolFromSmarts("CCCCC/C=C/C/C=C/C/C=C/C/C=C/C"),  # omega-6 (5,8,11,14)
        "epa": Chem.MolFromSmarts("CC/C=C/C/C=C/C/C=C/C/C=C/C/C=C/C"),  # omega-3 (5,8,11,14,17)
        "dha": Chem.MolFromSmarts("CC/C=C/C/C=C/C/C=C/C/C=C/C/C=C/C/C=C/C")  # omega-3 (4,7,10,13,16,19)
    }

    # Check for essential fatty acid patterns
    matched_patterns = []
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            matched_patterns.append(name)

    if not matched_patterns:
        # Check for methylene-interrupted double bond system
        methylene_interrupted = Chem.MolFromSmarts("C/C=C/CC/C=C/")
        if not mol.HasSubstructMatch(methylene_interrupted):
            return False, "Does not match essential fatty acid patterns"

    # Verify basic fatty acid properties
    if not is_phospholipid:
        formula = rdMolDescriptors.CalcMolFormula(mol)
        if not ('C' in formula and 'O2' in formula and 'H' in formula):
            return False, "Molecular formula not consistent with fatty acid"

    # Calculate molecular properties
    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rotatable_bonds < 3:
        return False, "Too rigid for an essential fatty acid"

    structure_type = "phospholipid" if is_phospholipid else "free fatty acid"
    pattern_desc = f" ({', '.join(matched_patterns)})" if matched_patterns else ""
    
    return True, f"Essential {structure_type} with {n_double_bonds} double bonds and {c_count} carbons{pattern_desc}"