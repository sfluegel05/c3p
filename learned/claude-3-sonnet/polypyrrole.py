"""
Classifies: CHEBI:38077 polypyrrole
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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

    # Kekulize the molecule to ensure proper aromatic perception
    Chem.Kekulize(mol, clearAromaticFlags=True)
    
    # Define improved pattern dictionary with more specific patterns
    pattern_dict = {
        # Basic aromatic pyrrole (excludes saturated rings)
        'pyrrole': '[nH,n]1:c:c:c:c1',
        
        # Metal-coordinated patterns (more specific)
        'metal_pyrrole': '[n]1(-[Fe,Mg,Zn,Co,Pd,B,#0]):c:c:c:c1',
        'porphyrin_metal': '[n]1:c:c:c2:c1:[n]:[#6]:c3:c2:[n]:[#6]:c:c3',
        
        # Porphyrin-type cores
        'porphyrin_core': '[nH,n]1:c:c:c2:c1:[nH,n]:c:c:c2',
        'chlorin_core': '[nH,n]1:c:c:c2:c1:[nH,n]:C:C:c2',
        
        # Direct connections between pyrroles
        'bipyrrole': '[nH,n]1:c:c:c:c1-[c]2:[nH,n]:c:c:c2',
        
        # Exclude these patterns (negative matches)
        'proline': '[N]1[CH2][CH2][CH2][CH]1',
        'indole': '[nH,n]1:c2:c:c:c:c:c2:c:c1'
    }
    
    # Convert patterns to RDKit molecules
    patterns = {}
    for name, smarts in pattern_dict.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            patterns[name] = patt
    
    # Count matches for each pattern
    matches = {}
    for name, patt in patterns.items():
        matches[name] = len(mol.GetSubstructMatches(patt))
    
    # Check for proline rings (exclude if only prolines found)
    if matches.get('proline', 0) > 0 and matches.get('pyrrole', 0) == 0:
        return False, "Contains only proline rings, not pyrroles"
    
    # Count total pyrrole units (including metal-coordinated)
    total_pyrroles = (
        matches.get('pyrrole', 0) + 
        matches.get('metal_pyrrole', 0) +
        matches.get('porphyrin_core', 0) * 2 +  # Count porphyrin cores as 2 pyrroles
        matches.get('chlorin_core', 0) * 2 +    # Count chlorin cores as 2 pyrroles
        matches.get('bipyrrole', 0) * 2         # Count bipyrroles as 2 pyrroles
    )
    
    # Special case: porphyrin with metal
    if matches.get('porphyrin_metal', 0) > 0:
        total_pyrroles = max(total_pyrroles, 4)  # Porphyrins always have 4 pyrroles
        
    # Check minimum requirement
    if total_pyrroles < 2:
        return False, f"Found only {total_pyrroles} pyrrole units, minimum of 2 required"
    
    # Determine structure type
    structure_types = []
    if matches.get('porphyrin_metal', 0) > 0 or matches.get('metal_pyrrole', 0) > 0:
        structure_types.append("metal-coordinated")
    if matches.get('porphyrin_core', 0) > 0:
        structure_types.append("porphyrin-like")
    if matches.get('bipyrrole', 0) > 0:
        structure_types.append("direct bipyrrole")
    if matches.get('chlorin_core', 0) > 0:
        structure_types.append("chlorin-type")
        
    # Additional validation
    ring_count = len(Chem.GetSymmSSSR(mol))
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    
    # Construct reason
    reason = f"Contains {total_pyrroles} pyrrole units"
    if structure_types:
        reason += f" ({', '.join(structure_types)})"
    if ring_count > 4 and aromatic_atoms > 12:
        reason += " in a complex polypyrrole system"
        
    return True, reason