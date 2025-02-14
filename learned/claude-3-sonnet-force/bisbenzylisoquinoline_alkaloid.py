"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    These compounds contain two benzylisoquinoline units linked by ether bridges or other connections.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic molecular properties check
    rings = rdMolDescriptors.CalcNumRings(mol)
    if rings < 6:
        return False, f"Too few rings ({rings}), need at least 6"

    # Count N and O atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 2:
        return False, f"Need at least 2 nitrogen atoms, found {n_count}"

    # More flexible pattern for isoquinoline/tetrahydroisoquinoline core
    # Allows for different oxidation states and substitution patterns
    isoquinoline_pattern = "[#7;R2]1[#6]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~[#6]~2~[#6]~1"
    isoquinoline_mol = Chem.MolFromSmarts(isoquinoline_pattern)
    
    if not isoquinoline_mol:
        return None, "Invalid SMARTS pattern"
        
    isoquinoline_matches = mol.GetSubstructMatches(isoquinoline_mol)
    num_isoquinoline = len(isoquinoline_matches)
    
    if num_isoquinoline < 2:
        return False, f"Found {num_isoquinoline} isoquinoline units, need at least 2"

    # Check for connecting groups between the isoquinoline units
    connecting_patterns = [
        "c-O-c",  # Ether bridge
        "c-O-C-O-c",  # Methylenedioxy bridge
        "c-c",  # Direct C-C connection
        "c-C-O-c"  # Extended ether bridge
    ]
    
    connection_count = 0
    for pattern in connecting_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            connection_count += len(mol.GetSubstructMatches(pat))

    if connection_count < 1:
        return False, "No suitable connecting bridges found"

    # Molecular weight check - adjusted range based on examples
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for this class"

    # Check for aromatic rings count
    aromatic_rings = len(mol.GetSubstructMatches(Chem.MolFromSmarts('a1aaaaa1')))
    if aromatic_rings < 4:
        return False, f"Too few aromatic rings ({aromatic_rings}), need at least 4"

    # Check for characteristic substituents
    substituents = {
        "methoxy": ("cOC", 0),
        "hydroxy": ("cO", 0),  # Simplified pattern
        "N-methyl or N-H": ("[#7;R2]", 0)  # Any ring nitrogen
    }
    
    for name, (pattern, _) in substituents.items():
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            substituents[name] = (pattern, len(mol.GetSubstructMatches(pat)))

    # Build detailed reason string
    substituent_str = ", ".join(f"{count} {name}" for name, (_, count) in substituents.items())
    reason = (f"Contains {num_isoquinoline} isoquinoline units with appropriate connecting bridges. "
             f"Found typical substituents ({substituent_str}). MW: {mol_wt:.1f}")
    
    return True, reason