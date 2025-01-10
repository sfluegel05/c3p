"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polypyrrole(smiles: str):
    """
    Determines if a molecule is a polypyrrole based on its SMILES string.
    A polypyrrole is a compound composed of two or more pyrrole units.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) indicating if molecule is a polypyrrole and reason
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define pyrrole patterns with their descriptions
    pattern_dict = {
        'basic_pyrrole': '[nH,n]1cccc1',  # Basic pyrrole ring
        'metal_pyrrole': '[n]1c(cc2)cc1',  # Metal-coordinated pyrrole
        'reduced_pyrrole': '[NH1,N]1CCCC1',  # Reduced pyrrole
        'fused_pyrrole': '[nH,n]1c2cccc1c2',  # Fused pyrrole system
        'porphyrin_core': '[n]1c(cc2)cc1[n]2',  # Porphyrin-like core
        'bipyrrole': '[nH,n]1cccc1-[c]1[c,n][c,n][c,n][c,n]1'  # Direct bipyrrole linkage
    }
    
    # Convert patterns to RDKit molecules and store valid ones
    valid_patterns = {}
    for name, pattern in pattern_dict.items():
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None:
            valid_patterns[name] = patt
            
    if not valid_patterns:
        return False, "Failed to create valid SMARTS patterns"
    
    # Find all pyrrole units
    total_matches = set()
    pattern_hits = {}
    
    for name, patt in valid_patterns.items():
        matches = mol.GetSubstructMatches(patt)
        pattern_hits[name] = len(matches)
        total_matches.update(matches)
    
    # Count unique pyrrole units
    num_pyrroles = len(total_matches)
    
    # Check for common metal centers in porphyrins and related structures
    metal_pattern = '[Fe,Mg,Zn,Co,Pd,B]'
    metal_mol = Chem.MolFromSmarts(metal_pattern)
    has_metal = False
    if metal_mol:
        has_metal = mol.HasSubstructMatch(metal_mol)
    
    # Analyze structure type
    structure_type = []
    if pattern_hits.get('porphyrin_core', 0) > 0:
        structure_type.append("porphyrin-like")
    if pattern_hits.get('bipyrrole', 0) > 0:
        structure_type.append("direct bipyrrole")
    if pattern_hits.get('fused_pyrrole', 0) > 0:
        structure_type.append("fused pyrrole system")
    if has_metal:
        structure_type.append("metal-coordinated")
        
    # Decision making
    if num_pyrroles < 2:
        return False, f"Only {num_pyrroles} pyrrole units found, minimum of 2 required"
        
    # Construct detailed reason
    reason = f"Contains {num_pyrroles} pyrrole units"
    if structure_type:
        reason += f" ({', '.join(structure_type)})"
        
    # Additional checks for specific cases
    ring_count = len(Chem.GetSymmSSSR(mol))
    if ring_count > 3 and num_pyrroles >= 4:
        reason += " in a complex polypyrrole system"
    
    return True, reason