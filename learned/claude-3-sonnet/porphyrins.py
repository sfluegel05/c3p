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

    # List of common metals found in porphyrins
    metals = {"Fe", "Mg", "Zn", "Co", "Ni", "Cu", "Pd", "Pt"}
    
    # Basic patterns for porphyrin core structures
    patterns = [
        # Basic porphyrin core (more flexible)
        "[#7]1~[#6]~[#6]~[#6]2~[#7]~[#6]~[#6]~[#6]3~[#7]~[#6]~[#6]~[#6]4~[#7]~[#6]~[#6]~[#6]~1~2~3~4",
        
        # Metalloporphyrin pattern
        "[#7]1~[#6]~[#6]~[#6]2~[#7]~[#6]~[#6]~[#6]3~[#7]~[#6]~[#6]~[#6]4~[#7]~[#6]~[#6]~[#6]~1~2~3~4~[Fe,Mg,Zn,Co,Ni,Cu,Pd,Pt]",
        
        # Pattern for reduced forms (chlorins/bacteriochlorins)
        "[#7]1~[#6]~[#6]~[#6]2~[#7]~[#6]~[#6]~[#6]3~[#7]~[#6]~[#6]~[#6]4~[#7]~[#6]~[#6]~[#6]~1~2~3~4",
        
        # Alternative pattern with explicit pyrrole rings
        "[nX2,nX3]1[cX3]~[cX3]~[#6]~[nX2,nX3]~[cX3]~[cX3]~[#6]~[nX2,nX3]~[cX3]~[cX3]~[#6]~[nX2,nX3]~[cX3]~[cX3]~[#6]~1"
    ]

    # Check for any matching pattern
    found_pattern = False
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt and mol.HasSubstructMatch(patt):
            found_pattern = True
            break
            
    if not found_pattern:
        return False, "No porphyrin core structure found"

    # Count nitrogens and check their environment
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if len(n_atoms) < 4:
        return False, "Insufficient nitrogen atoms for porphyrin structure"

    # Check for metal coordination
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() in metals]
    has_metal = len(metal_atoms) > 0

    # Ring analysis
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings found"

    # Find rings containing nitrogens
    n_containing_rings = 0
    large_ring_size = 0
    for ring in ring_info.AtomRings():
        ring_atoms = set(ring)
        n_atoms = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_atoms > 0:
            n_containing_rings += 1
            if len(ring) > large_ring_size:
                large_ring_size = len(ring)

    if n_containing_rings < 4:
        return False, "Insufficient number of nitrogen-containing rings"

    if large_ring_size < 16:
        return False, "Macrocyclic ring too small for porphyrin structure"

    # Count aromatic nitrogens vs total nitrogens
    aromatic_n = sum(1 for atom in n_atoms if atom.GetIsAromatic())
    
    # Determine the specific type of porphyrin
    if has_metal:
        metal = metal_atoms[0].GetSymbol()
        return True, f"Metalloporphyrin containing {metal}"
    elif aromatic_n >= 2:
        return True, "Free-base porphyrin"
    else:
        return True, "Reduced porphyrin (chlorin/bacteriochlorin)"