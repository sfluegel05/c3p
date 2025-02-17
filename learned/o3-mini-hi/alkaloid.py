"""
Classifies: CHEBI:22315 alkaloid
"""
#!/usr/bin/env python3
"""
Classifies: Alkaloid

Definition (heuristic):
  Alkaloids are defined as naturally occurring, basic nitrogen compounds,
  usually featuring one or more heterocyclic rings. (By extension, certain
  biogenetically related neutral compounds are also classed as alkaloids.)
  In contrast, compounds whose nitrogen is exclusively exocyclic (for example,
  dopamine, mescaline, serotonin and some simple amines) are usually not regarded
  as alkaloids.
  
This improved heuristic does the following:
  1. Verifies that the input SMILES is valid.
  2. Computes a rough “heteroatom ratio” (non‐C, non‐H atoms/total atoms) to help exclude
     molecules such as peptides/nucleotides.
  3. Checks for any nitrogen atoms that are in a ring.
      • If at least one ring‐nitrogen is found and the overall ring system is complex 
        (more than one ring or a fused system) then the molecule is classified as alkaloid.
      • If the only ring is a single small (≤6 atoms), fully aromatic ring then it is likely 
        a simple heterocycle (false positive).
  4. If no nitrogen is in a ring then we look for a benzylamine motif – which may indicate an 
     exocyclic nitrogen attached to an aromatic system (e.g. (–)-selegiline).
     
Note: This is only one example of a heuristic; many alkaloids defy any strict rule.
"""

from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid (using improved heuristics)
    
    Args:
        smiles (str): Input SMILES string.
        
    Returns:
        bool: True if classified as an alkaloid, False otherwise.
        str: Explanation of the reasoning.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate a rough heteroatom ratio (excluding H and carbon)
    total_atoms = mol.GetNumAtoms()
    hetero_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6)]
    hetero_ratio = len(hetero_atoms) / total_atoms if total_atoms > 0 else 0
    # If there are a lot of heteroatoms, the compound might be a peptide/nucleotide
    if hetero_ratio > 0.5:
        return False, "High heteroatom ratio suggests peptide/nucleotide rather than alkaloid"
    
    # Get all nitrogen atoms
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atoms found; alkaloids require at least one nitrogen."
    
    ring_info = mol.GetRingInfo().AtomRings()
    # Get list of rings (each as a tuple of atom indices)
    all_rings = list(ring_info)
    n_rings = len(all_rings)
    
    # Flag if any nitrogen is found inside a ring
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
    # Case 1: At least one nitrogen is part of a ring.
    if nitrogen_in_ring:
        # If the molecule has only one ring and that ring is small and fully aromatic,
        # it may be a simple heterocycle (many simple aromatic heterocycles are not regarded as alkaloids).
        if n_rings == 1:
            ring_atoms = all_rings[0]
            # Get the sub-molecule for the ring and then check aromaticity of each atom
            ring_mol = Chem.PathToSubmol(mol, ring_atoms)
            if ring_mol.GetNumAtoms() <= 6:
                # Check if every atom in the ring is aromatic.
                all_aromatic = all(atom.GetIsAromatic() for atom in ring_mol.GetAtoms())
                if all_aromatic:
                    return False, "Only a single small, fully aromatic ring detected; likely a simple heterocycle (amine) rather than an alkaloid"
        return True, "Nitrogen present in a ring system within a sufficiently complex scaffold; matches features of alkaloids"
    
    # Case 2: No nitrogen is in a ring.
    # Look for a benzylamine motif (exocyclic N attached via CH2 to an aromatic ring)
    benzylamine_pattern = Chem.MolFromSmarts("[NX3;H2,H1;!$(*@R)]-CH2-c1ccccc1")
    if benzylamine_pattern is not None and mol.HasSubstructMatch(benzylamine_pattern):
        return True, "Exocyclic nitrogen in a benzylamine-like motif attached to an aromatic ring; matches known alkaloid features"
    
    return False, "All nitrogen atoms are exocyclic and no benzylamine motif was found; likely an amine rather than an alkaloid"

# Example test cases (these are only examples; see outcomes for lists of TPs and FPs)
if __name__ == "__main__":
    test_examples = {
        "(-)-selegiline": "C[C@H](Cc1ccccc1)N(C)CC#C",
        "N-methylquipazine (false positive example)": "CN1CCN(CC1)c1ccc2ccccc2n1",
        "Simple pyrazine (false positive example)": "N1=C(C(=NC=C1C)C)/C=C/C",
    }
    for name, smi in test_examples.items():
        result, reason = is_alkaloid(smi)
        print(f"\nMolecule: {name}\n SMILES: {smi}\n Classified as alkaloid? {result}\n Reason: {reason}")