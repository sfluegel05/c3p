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

    # Check for the porphyrin macrocycle using a SMARTS pattern
    porphyrin_core = Chem.MolFromSmarts("[nH]1c(-[#6])c(-[#6])c2[nH]c(-[#6])c(-[#6])c3[nH]c(-[#6])c(-[#6])c4[nH]c(-[#6])c(-[#6])c1-2-3-4")
    if mol.HasSubstructMatch(porphyrin_core):
        return True, "Contains porphyrin macrocycle with four pyrrole rings connected by methine groups"

    # Alternative check for metal-coordinated variants (common in chlorophyll/heme)
    metal_porphyrin = Chem.MolFromSmarts("[n]1c(-[#6])c(-[#6])c2[n]c(-[#6])c(-[#6])c3[n]c(-[#6])c(-[#6])c4[n]c(-[#6])c(-[#6])c1-2-3-4")
    if mol.HasSubstructMatch(metal_porphyrin):
        return True, "Contains metalloporphyrin macrocycle structure"

    # Check for four pyrrole rings in a macrocycle
    pyrrole = Chem.MolFromSmarts("[nH]1cccc1")
    pyrrole_matches = mol.GetSubstructMatches(pyrrole)
    if len(pyrrole_matches) >= 4:
        # Check if all pyrroles are part of the same macrocycle
        ring_info = mol.GetRingInfo()
        macrocycle_found = False
        for ring in ring_info.AtomRings():
            if len(ring) >= 16:  # Approximate size for porphyrin macrocycle
                # Count pyrrole nitrogens in this ring
                nitro_in_ring = sum(1 for atom in ring if mol.GetAtomWithIdx(atom).GetAtomicNum() == 7)
                if nitro_in_ring >= 4:
                    macrocycle_found = True
                    break
        if macrocycle_found:
            return True, "Macrocycle with four pyrrole-like rings"

    return False, "Does not match porphyrin structural criteria"