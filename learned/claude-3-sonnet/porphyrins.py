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

    # Multiple SMARTS patterns to catch different porphyrin variants
    patterns = [
        # Basic porphyrin core (including reduced forms)
        "[nX2,nX3]1[cX3,CX4][cX3,CX4][cX3,CX4]2[cX3,CX4]3[nX2,nX3][cX3,CX4][cX3,CX4][cX3,CX4]4[cX3,CX4]5[nX2,nX3][cX3,CX4][cX3,CX4][cX3,CX4]6[cX3,CX4]7[nX2,nX3][cX3,CX4][cX3,CX4][cX3,CX4]8[cX3,CX4]1[cX3,CX4]2[cX3,CX4]3[cX3,CX4]4[cX3,CX4]5[cX3,CX4]6[cX3,CX4]7[cX3,CX4]8",
        # Metalloporphyrin pattern
        "[nX2]1(->[*]2<-[nX2]3[cX3,CX4][cX3,CX4][cX3,CX4]4[cX3,CX4]5[nX2](->[*]2)[cX3,CX4][cX3,CX4][cX3,CX4]6[cX3,CX4]7[nX2](->[*]2)[cX3,CX4][cX3,CX4][cX3,CX4]8[cX3,CX4]1[cX3,CX4]3[cX3,CX4]4[cX3,CX4]5[cX3,CX4]6[cX3,CX4]7[cX3,CX4]8)",
        # Chlorin/bacteriochlorin pattern (reduced porphyrins)
        "[nX2,nX3]1[cX3,CX4][cX3,CX4][cX3,CX4]2[cX3,CX4]3[nX2,nX3][CX4,cX3][CX4,cX3][cX3,CX4]4[cX3,CX4]5[nX2,nX3][cX3,CX4][cX3,CX4][cX3,CX4]6[cX3,CX4]7[nX2,nX3][CX4,cX3][CX4,cX3][cX3,CX4]8[cX3,CX4]1[cX3,CX4]2[cX3,CX4]3[cX3,CX4]4[cX3,CX4]5[cX3,CX4]6[cX3,CX4]7[cX3,CX4]8"
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

    # Count nitrogens - should have at least 4 for the core structure
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 4:
        return False, "Insufficient nitrogen atoms for porphyrin structure"

    # Analyze ring system
    ring_info = mol.GetRingInfo()
    ring_atoms = set()
    for ring in ring_info.AtomRings():
        ring_atoms.update(ring)
    
    if len(ring_atoms) < 20:  # Minimum atoms in ring system for porphyrin core
        return False, "Insufficient ring system size"

    # Check for metal coordination
    metals = ["Fe", "Mg", "Zn", "Co", "Ni", "Cu", "Pd", "Pt"]
    metal_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() in metals]
    has_metal = len(metal_atoms) > 0

    # Additional validation: check for connected ring system
    if not has_metal and n_count < 4:
        return False, "Free-base porphyrin requires at least 4 nitrogens"

    # Classify type of porphyrin
    if has_metal:
        metal = metal_atoms[0].GetSymbol()
        return True, f"Metalloporphyrin ({metal}) with tetrapyrrole macrocycle"
    else:
        return True, "Free-base porphyrin with tetrapyrrole macrocycle"