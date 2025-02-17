"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: proanthocyanidin
Definition: A flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans.
Heuristic:
  (1) The molecule must be valid.
  (2) It should have a molecular weight ≥480 Da.
  (3) It should contain at least 4 benzene rings (aromatic six‐membered rings).
  (4) It should have at least 6 hydroxyl (–OH) groups.
  (5) And it should either contain at least two hydroxyflavan-like substructures (using a relaxed catechin/epicatechin SMARTS)
      or show an inter-flavan linkage; the latter is tested by iteratively breaking candidate single bonds (outside rings)
      and verifying that each of the two resulting fragments has either (a) at least one benzene ring and a molecular weight 
      in the expected monomer range (~250–450 Da) OR (b) at least two benzene rings.
      
This revised code attempts to reduce both false positives (by not over accepting simple polyphenols)
and false negatives (by relaxing fragment criteria when the fragment weight is low).
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
      - And either contains at least two hydroxyflavan-like (catechin/epicatechin) substructures 
        (using a relaxed SMARTS that ignores chirality and allows flexibility)
        OR shows evidence of inter-flavan connectivity by bond fragmentation.
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule likely is a proanthocyanidin, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (1) Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 480:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da); expecting oligomers (≥480 Da)"
    
    # (2) Global aromatic ring count: Count aromatic six-membered (benzene-like) rings.
    ring_info = mol.GetRingInfo().AtomRings()
    benzene_rings = []
    for ring in ring_info:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            benzene_rings.append(set(ring))
    benzene_count = len(benzene_rings)
    if benzene_count < 4:
        return False, f"Not enough aromatic benzene rings (found {benzene_count}, need at least 4)"
    
    # (3) Count hydroxyl groups using a SMARTS for -OH (terminal oxygen with hydrogen).
    oh_smarts = "[OX2H]"
    oh_pattern = Chem.MolFromSmarts(oh_smarts)
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    n_oh = len(oh_matches)
    if n_oh < 6:
        return False, f"Not enough hydroxyl (-OH) groups (found {n_oh}, need at least 6)"
    
    # (4) First approach: Look for hydroxyflavan-like substructures.
    # We use a more relaxed SMARTS pattern that ignores chirality.
    # This pattern looks for a flavan core: an oxygen attached to a chain (without insisting on terminal -OH).
    hydroxyflavan_smarts = "O[C]1CCc2c(O)cc(O)c2C1"
    hydroxyflavan = Chem.MolFromSmarts(hydroxyflavan_smarts)
    # Using useChirality=False to allow variations in stereochemistry.
    flavan_matches = mol.GetSubstructMatches(hydroxyflavan, useChirality=False)
    if len(flavan_matches) >= 2:
        return True, ("Molecule contains at least two hydroxyflavan-like substructures based on a relaxed catechin/epicatechin pattern, "
                      "consistent with a proanthocyanidin (flavonoid oligomer).")
    
    # (5) Second approach: Look for inter-flavan connectivity by bond fragmentation.
    # For each bond that is a single bond and not in a ring, break it and evaluate the two fragments.
    inter_unit_link_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE or bond.IsInRing():
            continue
        # Try to fragment the molecule on this bond
        try:
            frag_mol = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)
        except Exception:
            continue
        frags = Chem.GetMolFrags(frag_mol, asMols=True)
        if len(frags) < 2:
            continue
        
        valid_frag_count = 0
        for frag in frags:
            frag_wt = rdMolDescriptors.CalcExactMolWt(frag)
            # Count aromatic benzene rings in the fragment.
            frag_ring_info = frag.GetRingInfo().AtomRings()
            frag_benzene_count = 0
            for ring in frag_ring_info:
                if len(ring) == 6 and all(frag.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    frag_benzene_count += 1
            # Define a fragment as a candidate flavan unit if either:
            # (a) it is in the expected monomer weight range (250-450 Da) and has at least 1 benzene ring,
            # or (b) it has at least 2 benzene rings.
            if (250 <= frag_wt <= 450 and frag_benzene_count >= 1) or (frag_benzene_count >= 2):
                valid_frag_count += 1
        if valid_frag_count >= 2:
            inter_unit_link_found = True
            break

    if not inter_unit_link_found:
        return False, ("No clear inter-flavan connectivity detected. Although overall weight, aromatic rings, and -OH counts pass, "
                       "the molecule lacks either two hydroxyflavan-like substructures or a candidate inter-flavan bond "
                       "whose fragmentation produces two plausible flavan units.")
    
    return True, ("Molecule has appropriate molecular weight, sufficient aromatic rings and hydroxyl groups, and either "
                  "contains multiple hydroxyflavan-like substructures or exhibits an inter-flavan linkage upon fragmentation; "
                  "consistent with a proanthocyanidin (flavonoid oligomer).")

# Example usage:
if __name__ == "__main__":
    # You can test with one of the known proanthocyanidins SMILES.
    test_smiles = ("COc1c(O)cc(cc1O)[C@H]1Oc2c(C[C@@H]1O)c(O)cc(O)c2"
                   "[C@H]1[C@H](O)[C@H](Oc2cc(O)cc(O)c12)c1cc(O)c(OC)c(O)c1")
    result, reason = is_proanthocyanidin(test_smiles)
    print(result, reason)