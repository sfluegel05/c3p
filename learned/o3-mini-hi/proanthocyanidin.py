"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: proanthocyanidin
Definition: A flavonoid oligomer obtained by the condensation of 
two or more units of hydroxyflavans.
Heuristic:
  (1) The molecule must be valid.
  (2) It should have a molecular weight ≥480 Da.
  (3) It should contain at least 4 benzene rings (aromatic six‐membered rings).
  (4) It should have at least 6 hydroxyl (–OH) groups.
  (5) And it should either contain at least two hydroxyflavan-like substructures (using a relaxed catechin/epicatechin SMARTS)
      or show an inter-flavan linkage; the latter is tested by breaking candidate single bonds (not in rings)
      and verifying that each fragment has at least 2 aromatic rings.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule might be a proanthocyanidin based on its SMILES string.
    Uses several heuristics:
      - Valid molecule with molecular weight ≥480 Da.
      - Contains at least 4 aromatic (benzene) rings.
      - Contains at least 6 hydroxyl (-OH) groups.
      - Either contains at least two hydroxyflavan-like (catechin/epicatechin) substructures 
        (using a relaxed SMARTS) OR shows evidence of inter-flavan connectivity by bond fragmentation.
        
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
        return False, f"Molecular weight too low ({mol_wt:.1f} Da); expecting oligomers (≥480 Da)"
    
    # Count aromatic six-membered (benzene-like) rings.
    ring_info = mol.GetRingInfo().AtomRings()
    benzene_rings = []
    for ring in ring_info:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            benzene_rings.append(set(ring))
    benzene_count = len(benzene_rings)
    if benzene_count < 4:
        return False, f"Not enough aromatic benzene rings (found {benzene_count}, need at least 4)"
    
    # Count hydroxyl groups using a SMARTS for -OH groups.
    oh_smarts = "[OX2H]"
    oh_pattern = Chem.MolFromSmarts(oh_smarts)
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    n_oh = len(oh_matches)
    if n_oh < 6:
        return False, f"Not enough hydroxyl (-OH) groups (found {n_oh}, need at least 6)"
    
    # First approach: Look for hydroxyflavan (catechin/epicatechin) core.
    # Use a slightly relaxed SMARTS that allows flexibility.
    # Original: "OC1CCc2c(O)cc(O)c2C1O"
    # We relax by not requiring the terminal -OH at one end (allowing for galloylation) 
    # and ignore chirality.
    hydroxyflavan_smarts = "OC1CCc2c(O)cc(O)c2C1"  # note: terminal -OH optional
    hydroxyflavan = Chem.MolFromSmarts(hydroxyflavan_smarts)
    flavan_matches = mol.GetSubstructMatches(hydroxyflavan, useChirality=False)
    if len(flavan_matches) >= 2:
        return True, ("Molecule contains at least two hydroxyflavan-like substructures based on the catechin/epicatechin pattern, "
                      "consistent with a proanthocyanidin (flavonoid oligomer).")
    
    # Second approach: Look for inter-flavan connectivity by bond fragmentation.
    # Look through all bonds that are single and not in a ring.
    inter_unit_link_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE or bond.IsInRing():
            continue
        # Use the bond index to fragment the molecule.
        try:
            frag_mol = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)
        except Exception:
            continue
        frags = Chem.GetMolFrags(frag_mol, asMols=True)
        # We require at least 2 fragments.
        if len(frags) < 2:
            continue
        # For each fragment, count its aromatic benzene rings.
        valid_frags = 0
        for frag in frags:
            ring_info_frag = frag.GetRingInfo().AtomRings()
            frag_benzene_count = 0
            for ring in ring_info_frag:
                if len(ring) == 6 and all(frag.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    frag_benzene_count += 1
            # We require a fragment to have at least 2 aromatic rings to be considered a flavan unit.
            if frag_benzene_count >= 2:
                valid_frags += 1
        if valid_frags >= 2:
            inter_unit_link_found = True
            break
    if not inter_unit_link_found:
        return False, ("No clear inter-flavan connectivity detected. Although weight, aromatic ring and -OH counts pass, "
                       "lack of at least two hydroxyflavan-like units or a bond whose fragmentation generates two fragments "
                       "with multiple aromatic rings suggests a single-unit flavonoid.")
    
    return True, ("Molecule has appropriate molecular weight, sufficient benzene rings and hydroxyl groups, and either "
                  "contains multiple hydroxyflavan-like substructures or exhibits an inter-flavan linkage when fragmented; "
                  "consistent with a proanthocyanidin (flavonoid oligomer).")

# Example usage (uncomment for testing):
# test_smiles = "COc1c(O)cc(cc1O)[C@H]1Oc2c(C[C@@H]1O)c(O)cc(O)c2[C@H]1[C@H](O)[C@H](Oc2cc(O)cc(O)c12)c1cc(O)c(OC)c(O)c1"
# result, reason = is_proanthocyanidin(test_smiles)
# print(result, reason)