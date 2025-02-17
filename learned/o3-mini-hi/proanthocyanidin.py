"""
Classifies: CHEBI:26267 proanthocyanidin
"""
"""
Classifies: proanthocyanidin
Definition: A flavonoid oligomer obtained by the condensation of two or more units of hydroxyflavans.
Heuristic (revised):
  (1) The molecule must be valid and have MW ≥500 Da.
  (2) It must contain at least 5 aromatic six‐membered (benzene) rings.
  (3) It must have at least 7 hydroxyl (-OH) groups.
  (4) It must either contain at least two flavan‐like substructures
      (using a relaxed SMARTS that targets a hydroxylated chroman/benzopyran core)
      OR show evidence for an inter‐flavan linkage. For the latter, only non‐ring, single bonds
      connecting sp³ carbons that are each attached to an aromatic neighbor are considered.
      Upon bond fragmentation, both fragments are required to have plausible monomer MWs (200–500 Da),
      at least one benzene ring, and a ring oxygen.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule might be a proanthocyanidin based on its SMILES string,
    using a set of heuristics to detect multiple flavan-like units.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a candidate proanthocyanidin, False otherwise.
        str: Explanation for decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # (1) Check molecular weight (set threshold to 500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da); expecting oligomers (≥500 Da)"
    
    # (2) Count aromatic benzene rings.
    rings = mol.GetRingInfo().AtomRings()
    benzene_rings = []
    for ring in rings:
        if len(ring) == 6 and all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            benzene_rings.append(set(ring))
    benzene_count = len(benzene_rings)
    if benzene_count < 5:
        return False, f"Not enough aromatic benzene rings (found {benzene_count}, need at least 5)"
    
    # (3) Count hydroxyl groups using SMARTS "[OX2H]".
    oh_smarts = "[OX2H]"
    oh_pat = Chem.MolFromSmarts(oh_smarts)
    oh_matches = mol.GetSubstructMatches(oh_pat)
    n_oh = len(oh_matches)
    if n_oh < 7:
        return False, f"Not enough hydroxyl (-OH) groups (found {n_oh}, need at least 7)"
    
    # (4a) Look for flavan-like substructures.
    # This SMARTS aims to catch a core of a flavan unit:
    # It requires a benzopyran (chroman) substructure with a fused benzene ring and at least one -OH on the B ring.
    # We ignore stereochemistry.
    flavan_smarts = "c1ccc(c(c1)O)C2COc3ccc(O)cc3C2"
    flavan_pat = Chem.MolFromSmarts(flavan_smarts)
    flavan_matches = mol.GetSubstructMatches(flavan_pat, useChirality=False)
    if len(flavan_matches) >= 2:
        return True, ("Molecule contains at least two flavan-like substructures based on a relaxed flavan SMARTS, "
                      "consistent with a proanthocyanidin (flavonoid oligomer).")
    
    # (4b) Second approach: Search for an inter-flavan bond.
    # Look for non-ring, single bonds connecting sp3 carbons that are each attached to at least one aromatic neighbor.
    interflavan_found = False
    for bond in mol.GetBonds():
        # Only consider single bonds not in rings.
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE or bond.IsInRing():
            continue
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Require both atoms are sp3 carbons.
        if atom1.GetAtomicNum() != 6 or atom2.GetAtomicNum() != 6:
            continue
        if atom1.GetHybridization().name != "SP3" or atom2.GetHybridization().name != "SP3":
            continue
        # Check that each atom has at least one neighbor that is aromatic.
        def has_aromatic_neighbor(atom):
            for nbr in atom.GetNeighbors():
                if nbr.GetIsAromatic():
                    return True
            return False
        if not has_aromatic_neighbor(atom1) or not has_aromatic_neighbor(atom2):
            continue

        # Try fragmenting the molecule by breaking this bond.
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
            # Only consider fragments in plausible monomer range.
            if not (200 <= frag_wt <= 500):
                continue
            # Count benzene rings in fragment.
            frag_rings = frag.GetRingInfo().AtomRings()
            frag_benzene = 0
            for ring in frag_rings:
                if len(ring) == 6 and all(frag.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    frag_benzene += 1
            # Check for at least one ring oxygen.
            ring_oxy = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 8 and atom.IsInRing())
            if frag_benzene >= 1 and ring_oxy >= 1:
                valid_fragments += 1
        if valid_fragments >= 2:
            interflavan_found = True
            break
            
    if interflavan_found:
        return True, ("Molecule meets molecular weight, aromatic ring, and -OH criteria and shows evidence for an inter-flavan linkage "
                      "by bond fragmentation; consistent with a proanthocyanidin (flavonoid oligomer).")
    
    return False, ("No clear flavan connectivity was detected: the molecule does not contain at least two identifiable flavan-like substructures, "
                   "nor is there an inter-flavan bond which, upon fragmentation, yields two candidate flavan monomers (200–500 Da, with aromatic ring and ring oxygen).")
    
# Example usage:
if __name__ == "__main__":
    # Example test SMILES from provided cases (one of the true positives)
    test_smiles = ("COc1c(O)cc(cc1O)[C@H]1Oc2c(C[C@@H]1O)c(O)cc(O)c2"
                   "[C@H]1[C@H](O)[C@H](Oc2cc(O)cc(O)c12)c1cc(O)c(OC)c(O)c1")
    result, explanation = is_proanthocyanidin(test_smiles)
    print(result, explanation)