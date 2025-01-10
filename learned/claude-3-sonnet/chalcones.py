"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: chalcones
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a 1,3-diphenylpropenone (benzylideneacetophenone) or its derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core chalcone patterns
    chalcone_patterns = [
        # Basic chalcone pattern (Ar-CH=CH-C(=O)-Ar)
        Chem.MolFromSmarts("[$([cR1,cR2]):1]!@[CH,C]=,:[CH,C]!@C(=O)!@[$([cR1,cR2]):2]"),
        
        # Dihydrochalcone pattern (Ar-CH2-CH2-C(=O)-Ar)
        Chem.MolFromSmarts("[$([cR1,cR2]):1]!@[CH2][CH2]!@C(=O)!@[$([cR1,cR2]):2]"),
        
        # Pattern for chalcones with modified double bond
        Chem.MolFromSmarts("[$([cR1,cR2]):1]!@[C,c]~[C,c]!@C(=O)!@[$([cR1,cR2]):2]")
    ]

    # Find matching pattern
    matching_pattern = None
    pattern_atoms = None
    
    for pattern in chalcone_patterns:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            matches = mol.GetSubstructMatches(pattern)
            if matches:
                pattern_atoms = matches[0]
                matching_pattern = pattern
                break
    
    if not pattern_atoms:
        return False, "No chalcone core structure found"

    # Basic validation
    ring_info = mol.GetRingInfo()
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 2:
        return False, "Insufficient number of rings"
    if ring_count > 8:
        return False, "Too many rings for typical chalcone"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, "Molecular weight outside typical chalcone range"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, "Too few carbons for chalcone structure"
    
    # Check for glycoside patterns (to filter out false positives)
    glycoside_pattern = Chem.MolFromSmarts("[OH1][CH1]1O[CH1][CH1][CH1][CH1]1")
    if mol.HasSubstructMatch(glycoside_pattern):
        sugar_count = len(mol.GetSubstructMatches(glycoside_pattern))
        if sugar_count > 1:
            return False, "Multiple sugar moieties present"

    # Determine chalcone type
    is_dihydro = mol.HasSubstructMatch(chalcone_patterns[1])
    chalcone_type = "dihydrochalcone" if is_dihydro else "chalcone"

    # Look for common substitution patterns
    substitution_patterns = {
        "hydroxy": Chem.MolFromSmarts("cO[H]"),
        "methoxy": Chem.MolFromSmarts("cOC"),
        "prenyl": Chem.MolFromSmarts("CC(C)=CCc"),
        "methylenedioxy": Chem.MolFromSmarts("OCOc")
    }
    
    substitutions = []
    for name, pattern in substitution_patterns.items():
        if pattern is not None and mol.HasSubstructMatch(pattern):
            substitutions.append(name)

    # Build classification reason
    if substitutions:
        substitution_desc = f" with {', '.join(substitutions)} groups"
        reason = f"Contains {chalcone_type} core structure{substitution_desc}"
    else:
        reason = f"Contains {chalcone_type} core structure"

    return True, reason