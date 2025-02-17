"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: proanthocyanidin
Definition: A flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans.
Heuristic:
  (1) The molecule must be valid.
  (2) It should have a molecular weight above ~480 Da (to capture slightly lower–weight borderline cases).
  (3) It should have at least 4 aromatic six-membered rings.
  (4) It should have several hydroxyl (–OH) groups (here, at least 6).
  (5) And ideally it should contain two or more flavan units as detected by a catechin–like substructure
      (SMARTS: OC1CCc2c(O)cc(O)c2C1O) OR show at least one bond connecting two different aromatic rings.
These heuristic filters are not perfect but are improved relative to the previous attempt.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule might be a proanthocyanidin based on its SMILES string.
    The heuristic:
      - must be valid,
      - has molecular weight ≥480 Da,
      - has at least 4 aromatic benzene rings (ring size=6),
      - has at least 6 hydroxyl groups,
      - and must either contain at least two flavan units (detected via a catechin core SMARTS)
        or show evidence of inter-flavan connectivity (at least one single bond linking two atoms
        that each belong to distinct aromatic six-membered rings).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule likely is a proanthocyanidin, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 480:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da); expecting oligomers (>480 Da)"
    
    # Count aromatic six-membered rings.
    # Get ring info and count rings of size 6 with all atoms aromatic.
    ring_info = mol.GetRingInfo().AtomRings()
    benzene_count = 0
    benzene_rings = []  # collect sets of atom indices for rings that are benzene-like
    for ring in ring_info:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            benzene_count += 1
            benzene_rings.append(set(ring))
    if benzene_count < 4:
        return False, f"Not enough aromatic benzene rings (found {benzene_count}, need at least 4)"
    
    # Count hydroxyl groups using a SMARTS for -OH.
    oh_smarts = "[OX2H]"
    oh_pattern = Chem.MolFromSmarts(oh_smarts)
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    n_oh = len(oh_matches)
    if n_oh < 6:
        return False, f"Not enough hydroxyl (-OH) groups (found {n_oh}, need at least 6)"
    
    # --- First approach: Look for a flavan (catechin/epicatechin) unit.
    # Many proanthocyanidins are built from catechin or epicatechin units.
    #
    # The following SMARTS tries to capture the catechin core:
    #    OC1CCc2c(O)cc(O)c2C1O
    # (Chirality is ignored so that both catechin and epicatechin match.)
    flavan_smarts = "OC1CCc2c(O)cc(O)c2C1O"
    flavan = Chem.MolFromSmarts(flavan_smarts)
    flavan_matches = mol.GetSubstructMatches(flavan, useChirality=False)
    # We require at least two nonoverlapping matches to indicate at least two flavan units.
    if len(flavan_matches) >= 2:
        return True, ("Molecule contains at least two catechin/epicatechin-like substructures "
                      "consistent with a proanthocyanidin (flavonoid oligomer).")
    
    # --- Second approach: Look for an inter-flavan bond.
    # We try to detect at least one single bond that connects atoms belonging
    # to two different aromatic rings (of set benzene_rings).
    inter_unit_link = False
    for bond in mol.GetBonds():
        # Check if the bond is a single bond and NOT in any ring (to avoid bonds within a fused system)
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE or bond.IsInRing():
            continue
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        # For each of the two atoms, get lists of benzene rings (by index in benzene_rings) that contain that atom.
        rings_a1 = [i for i, ring in enumerate(benzene_rings) if a1 in ring]
        rings_a2 = [i for i, ring in enumerate(benzene_rings) if a2 in ring]
        # Try to find a pair of rings (one from a1 and one from a2)
        # that are different (no shared atom)
        for i in rings_a1:
            for j in rings_a2:
                if i != j and benzene_rings[i].isdisjoint(benzene_rings[j]):
                    inter_unit_link = True
                    break
            if inter_unit_link:
                break
        if inter_unit_link:
            break
    if not inter_unit_link:
        return False, ("No inter-flavan linkage detected. Although weight, aromatic ring and -OH counts pass, "
                       "lack of a bond connecting two distinct aromatic units suggests a single-unit flavonoid.")
    
    return True, ("Molecule has appropriate molecular weight, multiple aromatic rings, several -OH groups, "
                  "and shows either two flavan units or an inter-flavan bond consistent with a proanthocyanidin (flavonoid oligomer).")
    
# Example usage (uncomment for testing):
# smiles = "COc1c(O)cc(cc1O)[C@H]1Oc2c(C[C@@H]1O)c(O)cc(O)c2[C@H]1[C@H](O)[C@H](Oc2cc(O)cc(O)c12)c1cc(O)c(OC)c(O)c1"
# result, reason = is_proanthocyanidin(smiles)
# print(result, reason)