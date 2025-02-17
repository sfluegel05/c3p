"""
Classifies: CHEBI:22315 alkaloid
"""
#!/usr/bin/env python3
"""
Classifies: Alkaloid

Definition (heuristic):
  Alkaloids are naturally occurring, basic nitrogen compounds,
  usually featuring one or more heterocyclic rings.
  (Some neutral compounds that are biogenetically related to basic alkaloids are also classed as alkaloids.)
  In contrast, compounds whose nitrogen is exclusively exocyclic (e.g. dopamine, mescaline, serotonin)
  are usually regarded as amines rather than alkaloids.
  
This implementation:
  1. Checks that the input SMILES is valid.
  2. Computes a rough “heteroatom ratio” to help filter out peptides/nucleotides.
  3. Searches for nitrogen atoms in rings. If at least one is found and the overall ring system is complex,
     the molecule is classified as an alkaloid. A single small fully aromatic ring is taken as a likely false positive.
  4. If no nitrogen is in a ring, it looks for a benzylamine motif.
"""

from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.

    Args:
        smiles (str): Input SMILES string.

    Returns:
        bool: True if the molecule is classified as an alkaloid, False otherwise.
        str: Explanation of the reasoning.
    """
    # Attempt to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute a rough heteroatom ratio (non-carbon, non-hydrogen atoms).
    total_atoms = mol.GetNumAtoms()
    hetero_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6)]
    hetero_ratio = len(hetero_atoms) / total_atoms if total_atoms > 0 else 0
    # If the heteroatom ratio is very high, the molecule might be a peptide or nucleotide.
    if hetero_ratio > 0.5:
        return False, "High heteroatom ratio suggests peptide/nucleotide rather than alkaloid"
    
    # Get all nitrogen atoms.
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found; alkaloids require at least one nitrogen."
    
    # Retrieve all rings from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    all_rings = list(ring_info)
    n_rings = len(all_rings)
    
    # Determine if at least one nitrogen is part of a ring.
    nitrogen_in_ring = False
    for atom in nitrogen_atoms:
        idx = atom.GetIdx()
        for ring in all_rings:
            if idx in ring:
                nitrogen_in_ring = True
                break
        if nitrogen_in_ring:
            break

    # Heuristic decision:
    # Case 1: At least one nitrogen is in a ring.
    if nitrogen_in_ring:
        # If the molecule has only one ring and that ring is small (≤6 atoms), check its aromaticity.
        if n_rings == 1:
            ring_atoms = all_rings[0]
            if len(ring_atoms) <= 6:
                # Instead of extracting a submol, check aromaticity directly from the parent molecule.
                all_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms)
                if all_aromatic:
                    return False, ("Only a single small, fully aromatic ring detected; "
                                   "likely a simple heterocycle (amine) rather than an alkaloid")
        return True, "Nitrogen present in a ring system within a sufficiently complex scaffold; matches features of alkaloids"
    
    # Case 2: No nitrogen is in a ring.
    # Check for a benzylamine-like motif: exocyclic N attached via CH2 to an aromatic ring.
    benzylamine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(*@R)]-CH2-c1ccccc1")
    if benzylamine_pattern is not None and mol.HasSubstructMatch(benzylamine_pattern):
        return True, ("Exocyclic nitrogen in a benzylamine-like motif attached to an aromatic ring; "
                      "matches known alkaloid features")
    
    return False, "All nitrogen atoms are exocyclic and no benzylamine motif was found; likely an amine rather than an alkaloid"

# Example test cases (these serve to illustrate the results)
if __name__ == "__main__":
    test_examples = {
        "(-)-selegiline": "C[C@H](Cc1ccccc1)N(C)CC#C",
        "Example with chloride ion": "[Cl-].CCn1cc[n+](C)c1",
        "Simple pyrazine (false positive example)": "N1=C(C(=NC=C1C)C)/C=C/C",
    }
    for name, smi in test_examples.items():
        result, reason = is_alkaloid(smi)
        print(f"\nMolecule: {name}\n SMILES: {smi}\n Classified as alkaloid? {result}\n Reason: {reason}")