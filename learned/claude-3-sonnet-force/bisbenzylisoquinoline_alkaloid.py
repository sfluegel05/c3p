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
    These compounds contain two benzylisoquinoline units linked by ether bridges.

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
    if o_count < 1:
        return False, f"Need oxygen atoms for ether bridges, found {o_count}"

    # Pattern for tetrahydroisoquinoline core (more flexible)
    thiq_pattern = """
        [#6]~1~[#6]~[#6]~[#7X3]~[#6]~[#6]~2~[#6]~[#6]~[#6]~[#6]~[#6]~12
    """
    thiq_mol = Chem.MolFromSmarts(thiq_pattern)
    
    if thiq_mol is None:
        return False, "Invalid SMARTS pattern"
        
    thiq_matches = mol.GetSubstructMatches(thiq_mol)
    num_thiq = len(thiq_matches)
    
    if num_thiq < 2:
        return False, f"Found {num_thiq} tetrahydroisoquinoline units, need at least 2"

    # Look for ether bridges between aromatic rings
    ether_patterns = [
        "c-O-c",  # Direct ether bridge
        "c-O-C-c", # Extended ether bridge
    ]
    
    ether_count = 0
    for pattern in ether_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            ether_count += len(mol.GetSubstructMatches(pat))

    if ether_count < 1:
        return False, "No suitable ether bridges found"

    # Check for characteristic substituents
    substituents = {
        "methoxy": ("cOC", 0),
        "hydroxy": ("cO[H]", 0),
        "N-methyl": ("[NX3]-C", 0)
    }
    
    for name, (pattern, _) in substituents.items():
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            substituents[name] = (pattern, len(mol.GetSubstructMatches(pat)))
    
    # Molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 450 or mol_wt > 800:
        return False, f"Molecular weight {mol_wt:.1f} outside typical range for this class"

    # Check for benzyl groups
    benzyl_pattern = Chem.MolFromSmarts("c1ccccc1-[CH2]")
    if benzyl_pattern:
        benzyl_count = len(mol.GetSubstructMatches(benzyl_pattern))
        if benzyl_count < 2:
            return False, f"Found only {benzyl_count} benzyl groups, need at least 2"

    # Build detailed reason string
    substituent_str = ", ".join(f"{count} {name}" for name, (_, count) in substituents.items())
    reason = (f"Contains {num_thiq} tetrahydroisoquinoline units, {ether_count} ether bridges, "
             f"and typical substituents ({substituent_str}). MW: {mol_wt:.1f}")
    
    return True, reason