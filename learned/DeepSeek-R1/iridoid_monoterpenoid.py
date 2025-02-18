"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
"""
Classifies: iridoid monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.
    Criteria: Cyclopentane ring fused to a six-membered oxygen-containing ring.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Get all rings in the molecule
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()
    
    # Find cyclopentane (5-membered all carbon)
    cyclo_rings = []
    for ring in rings:
        if len(ring) == 5:
            all_carbon = all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in ring)
            if all_carbon:
                cyclo_rings.append(ring)
    
    # Find six-membered rings with oxygen
    oxy_rings = []
    for ring in rings:
        if len(ring) == 6:
            has_oxygen = any(mol.GetAtomWithIdx(a).GetAtomicNum() == 8 for a in ring)
            if has_oxygen:
                oxy_rings.append(ring)
    
    # Check for fused cyclopentane and oxygen ring
    fused = False
    for cyclo in cyclo_rings:
        for oxy in oxy_rings:
            shared = set(cyclo) & set(oxy)
            if len(shared) >= 2:
                # Check if any two shared atoms are consecutive in either ring
                # Check cyclo ring
                for i in range(len(cyclo)):
                    a, b = cyclo[i], cyclo[(i+1)%5]
                    if {a, b}.issubset(shared):
                        fused = True
                        break
                if not fused:
                    # Check oxy ring
                    for i in range(len(oxy)):
                        a, b = oxy[i], oxy[(i+1)%6]
                        if {a, b}.issubset(shared):
                            fused = True
                            break
                if fused:
                    break
            if fused:
                break
        if fused:
            break
    
    if not fused:
        return False, "No fused cyclopentane and six-membered oxygen ring"
    
    # Optional: Check for isoprene-derived structure (monoterpenoid)
    # This is approximate; monoterpenoids have ~10 carbons but substituents vary
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:  # Adjust based on examples
        return False, "Too few carbons for monoterpenoid"
    
    return True, "Fused cyclopentane and six-membered oxygen-containing ring"