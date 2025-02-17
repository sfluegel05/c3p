"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: proanthocyanidin
Definition: A flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans.
Heuristic (revised):
  (1) The molecule must be valid and have MW ≥480 Da.
  (2) It must contain at least 4 aromatic six‐membered (benzene) rings.
  (3) It must have at least 6 hydroxyl (-OH) groups.
  (4) And it must either contain at least two flavan-like substructures
      (using a SMARTS that requires a hydroxylated chroman ring with at least one ring oxygen)
      OR show evidence for an inter-flavan linkage, as determined by bond fragmentation where
      each fragment (candidate flavan monomer) has a plausible monomer MW (250–450 Da),
      at least one benzene ring, and at least one oxygen atom in a ring.
      
This revised heuristic aims to reduce false positives from simple polyphenols and to catch cases
where interflavan connectivity is evident.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule might be a proanthocyanidin based on its SMILES string.
    
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

    # (1) Check molecular weight threshold
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 480:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da); expecting oligomers (≥480 Da)"

    # (2) Count aromatic benzene rings.
    rings = mol.GetRingInfo().AtomRings()
    benzene_rings = []
    for ring in rings:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            benzene_rings.append(set(ring))
    benzene_count = len(benzene_rings)
    if benzene_count < 4:
        return False, f"Not enough aromatic benzene rings (found {benzene_count}, need at least 4)"

    # (3) Count hydroxyl groups using SMARTS "[OX2H]"
    oh_smarts = "[OX2H]"
    oh_pat = Chem.MolFromSmarts(oh_smarts)
    oh_matches = mol.GetSubstructMatches(oh_pat)
    n_oh = len(oh_matches)
    if n_oh < 6:
        return False, f"Not enough hydroxyl (-OH) groups (found {n_oh}, need at least 6)"

    # (4) Look for flavan-like substructures.
    # This SMARTS is designed to capture the core of a flavan unit.
    # It requires a saturated heterocycle (chroman) with one ether oxygen and an attached aromatic ring (B-ring)
    # substituted with at least one hydroxyl.
    flavan_smarts = "O[C@H]1CC[C@H](O)[C@@H]1c1ccc(O)cc1"  
    # We ignore stereochemistry to allow flexibility.
    flavan_pat = Chem.MolFromSmarts(flavan_smarts)
    flavan_matches = mol.GetSubstructMatches(flavan_pat, useChirality=False)
    if len(flavan_matches) >= 2:
        return True, ("Molecule contains at least two flavan-like substructures based on a specialized flavan SMARTS, "
                      "consistent with a proanthocyanidin (flavonoid oligomer).")
    
    # (5) Second approach: Evaluate inter-flavan connectivity by bond fragmentation.
    # For each candidate single non-ring bond, break and check if both fragments are plausible flavan monomers.
    interflavan_found = False
    for bond in mol.GetBonds():
        # Only consider single bonds not in rings.
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE or bond.IsInRing():
            continue
        try:
            frag_mol = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)
        except Exception:
            continue
        frags = Chem.GetMolFrags(frag_mol, asMols=True)
        if len(frags) < 2:
            continue
        
        valid_fragments = 0
        for frag in frags:
            frag_wt = rdMolDescriptors.CalcExactMolWt(frag)
            # Count benzene rings in fragment
            frag_rings = frag.GetRingInfo().AtomRings()
            frag_benzene = 0
            for r in frag_rings:
                if len(r) == 6 and all(frag.GetAtomWithIdx(idx).GetIsAromatic() for idx in r):
                    frag_benzene += 1
            # Also check for presence of at least one oxygen in a ring (a proxy for the chroman oxygen)
            ring_oxy = 0
            for atom in frag.GetAtoms():
                if atom.GetAtomicNum() == 8 and atom.IsInRing():
                    ring_oxy += 1

            # Define a fragment as a candidate flavan unit if:
            # (a) Its weight is between 250 and 450 Da,
            # (b) It has at least one aromatic (benzene) ring,
            # (c) And it has at least one ring oxygen.
            if 250 <= frag_wt <= 450 and frag_benzene >= 1 and ring_oxy >= 1:
                valid_fragments += 1

        if valid_fragments >= 2:
            interflavan_found = True
            break

    if not interflavan_found:
        return False, ("No clear inter-flavan connectivity detected. Although overall weight, aromatic rings, and -OH counts pass, "
                       "the molecule lacks two recognizable flavan-like substructures or a candidate inter-flavan bond whose "
                       "fragmentation yields plausible flavan units (each with appropriate weight, benzene ring(s), and ring oxygen).")
    
    return True, ("Molecule meets molecular weight, aromatic ring, and -OH criteria and either contains multiple flavan-like substructures "
                  "or shows a plausible inter-flavan connectivity upon fragmentation; consistent with a proanthocyanidin (flavonoid oligomer).")

# Example usage:
if __name__ == "__main__":
    # Test with one known proanthocyanidin SMILES (taken from one of the examples)
    test_smiles = ("COc1c(O)cc(cc1O)[C@H]1Oc2c(C[C@@H]1O)c(O)cc(O)c2"
                   "[C@H]1[C@H](O)[C@H](Oc2cc(O)cc(O)c12)c1cc(O)c(OC)c(O)c1")
    result, explanation = is_proanthocyanidin(test_smiles)
    print(result, explanation)