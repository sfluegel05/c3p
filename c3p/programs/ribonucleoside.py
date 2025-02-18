"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: any nucleoside where the sugar component is D-ribose (ribonucleoside)

Definition:
A ribonucleoside is a nucleoside that has its nucleobase linked via an N-glycosidic bond to a D-ribofuranose unit.
The ribofuranose is a five-membered ring (4 carbons and 1 oxygen) and possesses an exocyclic CH2OH group 
(at the 5'-position) attached to one of its carbons. Additionally, one ring carbon is directly bonded 
to an aromatic nitrogen (representing the nucleobase linkage).

Improvements over previous attempts:
  1. Reject compounds containing phosphorus (avoid nucleotides).
  2. Look for five-membered rings having exactly 1 oxygen and 4 carbons.
  3. For each candidate ring, check that one of its carbons is attached to an exocyclic CH2OH group.
  4. Also require that at least one ring carbon is attached (outside the ring) to an aromatic nitrogen.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside (nucleoside with D-ribose sugar)
    by checking for:
      - absence of phosphorus,
      - presence of a five-membered ring (4 carbons, 1 oxygen),
      - at least one ring carbon bearing an exocyclic CH2OH group,
      - and a linkage (via an aromatic nitrogen) to a nucleobase.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a ribonucleoside, else False.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Reject molecules containing phosphorus (e.g. nucleotides)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; likely a nucleotide rather than a nucleoside"

    # For proper substituent perception, add explicit hydrogens.
    molH = Chem.AddHs(mol)
    ring_info = molH.GetRingInfo()

    # Helper function: determine if an atom is a CH2 (methylene) that is part of a CH2OH group.
    def is_CH2OH(candidate, attached_from_ring):
        # Candidate must be a carbon (atomic number 6)
        if candidate.GetAtomicNum() != 6:
            return False
        # Check that candidate is not part of the sugar ring (exocyclic)
        # Count explicit hydrogens (if not provided, GetTotalNumHs() can be used)
        # We expect a CH2 to have two hydrogens.
        if candidate.GetTotalNumHs() != 2:
            return False
        # It should be attached to exactly one oxygen that is itself an -OH group.
        oxygen_neighbors = [nbr for nbr in candidate.GetNeighbors() 
                            if nbr.GetAtomicNum() == 8 and nbr.GetIdx() != attached_from_ring]
        if len(oxygen_neighbors) != 1:
            return False
        # Check that this oxygen has at least one hydrogen (i.e. is -OH)
        oxygen = oxygen_neighbors[0]
        if oxygen.GetTotalNumHs() < 1:
            return False
        return True

    # Iterate over each ring in the molecule.
    for ring in ring_info.AtomRings():
        # We only consider rings with exactly five atoms.
        if len(ring) != 5:
            continue
        
        # Gather atoms of the ring (by indices in the current molH)
        ring_atoms = [molH.GetAtomWithIdx(idx) for idx in ring]
        # Count how many are oxygen (atomic number 8) and carbon (atomic number 6)
        num_oxygens = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        num_carbons = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
        if num_oxygens != 1 or num_carbons != 4:
            continue  # This ring does not match the ribofuranose pattern.

        # Now, check for the exocyclic CH2OH group.
        ch2oh_found = False
        # Also check for nucleobase linkage via aromatic nitrogen.
        nucleobase_attached = False
        
        # For every atom in the ring:
        for idx in ring:
            ring_atom = molH.GetAtomWithIdx(idx)
            # Check neighbors not in the ring.
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # Skip atoms inside the ring.
                # Check if this neighbor is a CH2OH group.
                # We require that the exocyclic candidate attached to the ring is a carbon that qualifies as CH2OH.
                if nbr.GetAtomicNum() == 6 and is_CH2OH(nbr, ring_atom.GetIdx()):
                    ch2oh_found = True
                # Also check if a neighbor is an aromatic nitrogen; this can mark nucleobase attachment.
                if nbr.GetAtomicNum() == 7 and nbr.GetIsAromatic():
                    nucleobase_attached = True
            # If we already found a nucleobase link on one of the ring atoms,
            # we can also check for the CH2OH on some (possibly different) ring atom.
            # (Do not short-circuit until both are true.)
        
        # If both features are found in this candidate ribofuranose ring, we classify as ribonucleoside.
        if ch2oh_found and nucleobase_attached:
            return True, ("Found five-membered ribofuranose ring (1 O and 4 C) with an exocyclic CH2OH substituent "
                          "and a ring carbon attached to an aromatic nitrogen (nucleobase linkage)")
    
    return False, "No ribose moiety (with CH2OH and nucleobase linkage) found"


# Example usage and testing (if run as __main__)
if __name__ == "__main__":
    examples = {
        "1-methyladenosine": "Cn1cnc2n(cnc2c1=N)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O",
        "3,4-dihydrozebularine": "OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)N1C=CCNC1=O",
        "nucleocidin": "NC1=C2N=CN([C@@H]3O[C@](F)(COS(N)(=O)=O)[C@@H](O)[C@H]3O)C2=NC=N1",
        "cytidine": "Nc1ccn([C@@H]2O[C@H](CO)[C@@H](O)[C@H]2O)c(=O)n1",
        "aminodeoxyfutalosine": "Nc1ncnc2n(cnc12)[C@@H]1O[C@H](CCC(=O)c2cccc(c2)C(O)=O)[C@@H](O)[C@H]1O",
        "5'-deoxyadenosine (should fail ribonucleoside)": "C[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(N)ncnc12",
        "Adenosylcobinamide phosphate (should be rejected)": "C[C@H](CNC(=O)CC[C@]1(C)[C@@H](CC(N)=O)[C@H]2N3C1=C(C)C1=[N+]4C(=CC5=[N+]6C(=C(C)C7=[N+]([C@]2(C)[C@@](C)(CC(N)=O)[C@@H]7CCC(N)=O)[Co--]346C[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@](C)(CC(N)=O)[C@@H]5CCC(N)=O)C(C)(C)[C@@H]1CCC(N)=O)OP(O)(O)=O"
    }
    
    for name, smi in examples.items():
        res, reason = is_ribonucleoside(smi)
        print(f"{name}: {res} ({reason})")