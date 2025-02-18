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

    # Improved SMARTS pattern for porphyrin core (four pyrroles connected by methine bridges)
    porphyrin_core = Chem.MolFromSmarts("[nH]1c2c([nH]c3c([nH]c4c([nH]c1cc4)cc3)cc2)cc1")
    if mol.HasSubstructMatch(porphyrin_core):
        return True, "Contains porphyrin macrocycle with four pyrrole rings connected by methine groups"

    # Check for metalloporphyrins (common in heme/chlorophyll)
    metal_porphyrin = Chem.MolFromSmarts("[Mg,Fe,Zn,Co,Pd,Pt]@[n]1c2c([n]c3c([n]c4c([n]c1cc4)cc3)cc2)cc1")
    if mol.HasSubstructMatch(metal_porphyrin):
        return True, "Contains metalloporphyrin core structure"

    # Check for macrocycle with 16 atoms (porphine skeleton size)
    ring_info = mol.GetRingInfo()
    macrocycle = False
    for ring in ring_info.AtomRings():
        if len(ring) == 16:
            # Check for four nitrogens in the macrocycle
            nitro_count = sum(1 for atom in ring if mol.GetAtomWithIdx(atom).GetAtomicNum() == 7)
            if nitro_count == 4:
                # Verify each nitrogen is in a pyrrole-like 5-membered ring
                pyrrole_nitro = 0
                for atom in ring:
                    a = mol.GetAtomWithIdx(atom)
                    if a.GetAtomicNum() == 7:
                        # Check if nitrogen is in a 5-membered aromatic ring
                        for bond in a.GetBonds():
                            if bond.GetBondType() == Chem.BondType.AROMATIC:
                                neigh = bond.GetOtherAtomIdx(atom)
                                # Check if the ring is 5-membered
                                for r in ring_info.AtomRings():
                                    if atom in r and neigh in r and len(r) == 5:
                                        pyrrole_nitro +=1
                                        break
                if pyrrole_nitro >=4:
                    macrocycle = True
                    break
    if macrocycle:
        return True, "16-membered macrocycle with four pyrrole-like nitrogens"

    return False, "Does not match porphyrin structural criteria"