"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:26355 porphyrins
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    Porphyrins have a macrocyclic structure with four pyrrole nuclei connected by methine groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Validate and create SMARTS patterns with error handling
    try:
        # Core porphyrin pattern (16-membered macrocycle with four pyrrole subunits)
        porphyrin_core = Chem.MolFromSmarts("c1c2nc3c4nc1cccc1c3nc2c(c1)c4")
        # Metalloporphyrin pattern (metal coordinated to four nitrogens)
        metal_porphyrin = Chem.MolFromSmarts("[Mg,Fe,Zn,Co,Pd,Pt]~[n]~[n]~[n]~[n]")
    except:
        return False, "SMARTS pattern error"

    if porphyrin_core is None or metal_porphyrin is None:
        return False, "SMARTS initialization failed"

    # Check core porphyrin structure
    if mol.HasSubstructMatch(porphyrin_core):
        return True, "Contains porphyrin macrocycle core"

    # Check for metalloporphyrins (common in biological systems)
    if mol.HasSubstructMatch(metal_porphyrin):
        # Verify the metal is part of a macrocycle
        ring_info = mol.GetRingInfo()
        for ring in ring_info.AtomRings():
            if len(ring) >= 16:  # Core porphyrin ring size
                metal_in_ring = any(mol.GetAtomWithIdx(a).GetAtomicNum() in {12, 26, 30, 27, 46, 78} for a in ring)
                if metal_in_ring:
                    return True, "Metalloporphyrin structure detected"

    # Macrocycle verification fallback
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if 16 <= len(ring) <= 24:  # Allow some substitution variations
            nitro_count = sum(1 for a in ring if mol.GetAtomWithIdx(a).GetAtomicNum() == 7)
            if nitro_count == 4:
                # Verify each nitrogen is in a pyrrole-like 5-membered ring
                pyrrole_nitro = 0
                for atom in ring:
                    a = mol.GetAtomWithIdx(atom)
                    if a.GetAtomicNum() == 7:
                        # Check for participation in a 5-membered ring
                        for r in ring_info.AtomRings():
                            if atom in r and len(r) == 5:
                                pyrrole_nitro += 1
                                break
                if pyrrole_nitro >= 4:
                    return True, "Macrocycle with four pyrrole-like nitrogens"

    return False, "Does not match porphyrin structural criteria"