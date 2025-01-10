"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: porphyrins
Natural pigments containing a fundamental skeleton of four pyrrole nuclei united through 
the alpha-positions by four methine groups to form a macrocyclic structure.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_porphyrin, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic porphyrin patterns - more flexible than before
    patterns = [
        # Basic porphyrin core with variations in bond types
        "[nX2,nX3]1[cX3,CX4]~[cX3,CX4]~[cX3,CX4]~[#6]~[nX2,nX3]~[cX3,CX4]~[cX3,CX4]~[cX3,CX4]~[#6]~[nX2,nX3]~[cX3,CX4]~[cX3,CX4]~[cX3,CX4]~[#6]~[nX2,nX3]~[cX3,CX4]~[cX3,CX4]~[cX3,CX4]~1",
        # Alternative pattern for metalloporphyrins
        "[#7]1~[#6]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~1",
        # Reduced forms (chlorins/bacteriochlorins)
        "[nX2,nX3]1[cX3,CX4]~[cX3,CX4]~[cX3,CX4]~[#6]~[nX2,nX3]~[CX4]~[CX4]~[cX3,CX4]~[#6]~[nX2,nX3]~[cX3,CX4]~[cX3,CX4]~[cX3,CX4]~[#6]~[nX2,nX3]~[cX3,CX4]~[cX3,CX4]~[cX3,CX4]~1"
    ]

    # Check for matching patterns
    found_pattern = False
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            found_pattern = True
            break
            
    if not found_pattern:
        return False, "No porphyrin core structure found"

    # Count nitrogens
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 4:
        return False, "Insufficient nitrogen atoms for porphyrin structure"

    # Ring analysis
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings found"

    # Find the largest ring system that contains nitrogens
    max_ring_size = 0
    n_containing_rings = 0
    for ring in ring_info.AtomRings():
        ring_atoms = set(ring)
        n_atoms = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_atoms > 0:
            n_containing_rings += 1
            max_ring_size = max(max_ring_size, len(ring))

    if n_containing_rings < 4:
        return False, "Insufficient number of nitrogen-containing rings"
    
    if max_ring_size < 16:  # Minimum size for the macrocycle
        return False, "Ring system too small for porphyrin"

    # Check for metals (optional feature)
    metals = {"Fe", "Mg", "Zn", "Co", "Ni", "Cu", "Pd", "Pt"}
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() in metals]
    has_metal = len(metal_atoms) > 0

    # Additional structure validation
    aromatic_n = sum(1 for atom in mol.GetAtoms() 
                    if atom.GetAtomicNum() == 7 and atom.GetIsAromatic())
    saturated_n = n_count - aromatic_n

    if has_metal:
        metal = metal_atoms[0].GetSymbol()
        return True, f"Metalloporphyrin ({metal}) with tetrapyrrole macrocycle"
    elif aromatic_n >= 2:
        return True, "Free-base porphyrin with tetrapyrrole macrocycle"
    elif saturated_n >= 4:
        return True, "Reduced porphyrin (chlorin/bacteriochlorin) structure"
    else:
        return True, "Porphyrin-like structure with tetrapyrrole core"