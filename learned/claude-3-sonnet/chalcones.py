"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: chalcones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # More specific chalcone core patterns
    chalcone_patterns = [
        # Classic chalcone pattern with explicit connection points
        Chem.MolFromSmarts("[$(c:c):1]!@[CH]=[CH]!@C(=O)!@[$(c:c):2]"), 
        
        # Dihydrochalcone pattern with explicit connection points
        Chem.MolFromSmarts("[$(c:c):1]!@[CH2][CH2]!@C(=O)!@[$(c:c):2]"),
        
        # Pattern for substituted chalcones
        Chem.MolFromSmarts("[$(c:c):1]!@[C]=[C]!@C(=O)!@[$(c:c):2]"),
        
        # Pattern for chalcones with fused rings
        Chem.MolFromSmarts("[$(c:c:c):1]!@[CH]=[CH]!@C(=O)!@[$(c:c):2]")
    ]

    # Find matching pattern and get mapped atoms
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

    # Verify aromatic rings are properly connected to the core
    ring_info = mol.GetRingInfo()
    aromatic_rings = []
    
    for ring in ring_info.AtomRings():
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings.append(ring)

    if len(aromatic_rings) < 2:
        return False, "Insufficient number of aromatic rings"

    # Check for required substitution patterns
    substitution_patterns = {
        "hydroxy": Chem.MolFromSmarts("cO"),
        "methoxy": Chem.MolFromSmarts("cOC"),
        "prenyl": Chem.MolFromSmarts("CC(C)=CCc"),
        "methylenedioxy": Chem.MolFromSmarts("OCOc")
    }
    
    substitutions = []
    for name, pattern in substitution_patterns.items():
        if pattern is not None and mol.HasSubstructMatch(pattern):
            substitutions.append(name)

    if not substitutions:
        return False, "No typical chalcone substitution patterns found"

    # Additional validation
    ketone_pattern = Chem.MolFromSmarts("C(=O)")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, "Too few carbons for chalcone structure"
    
    if o_count < 1:
        return False, "No oxygen atoms found"

    # Determine chalcone type
    is_dihydro = mol.HasSubstructMatch(chalcone_patterns[1])
    chalcone_type = "dihydrochalcone" if is_dihydro else "chalcone"
    
    # Build classification reason
    substitution_desc = f" with {', '.join(substitutions)} groups"
    reason = f"Contains {chalcone_type} core structure{substitution_desc}"

    return True, reason