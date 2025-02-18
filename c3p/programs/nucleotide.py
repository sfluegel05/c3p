"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: Nucleotide
A nucleotide is defined as a nucleoside phosphate resulting from the condensation 
of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.
This function verifies:
   - Presence of a furanose sugar ring (5-membered ring with 1 oxygen and 4 carbons)
     and considers exocyclic atoms attached to the sugar ring (e.g. the 5'-CH2 group).
   - A phosphate (P) is connected via an oxygen (O) to an atom in the sugar region.
   - A nucleobase-like aromatic heterocycle (aromatic ring with ≥ 2 nitrogen atoms)
     that is directly connected to the sugar region.
"""
from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    
    A nucleotide should contain a nucleoside portion (a sugar moiety attached to a 
    nucleobase) with a phosphate group esterified at a sugar hydroxy position.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a nucleotide, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    sugar_ring_indices = None
    # Step 1: Identify a furanose sugar ring: a 5-membered ring with exactly 1 oxygen and 4 carbons.
    for ring in ring_info.AtomRings():
        if len(ring) == 5:  # candidate for furanose ring
            symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
            if symbols.count('O') == 1 and symbols.count('C') == 4:
                sugar_ring_indices = set(ring)
                break
    if sugar_ring_indices is None:
        return False, "No furanose sugar ring (5-membered ring with 1 oxygen and 4 carbons) found"
    
    # Expand sugar region: include atoms not in the ring that are directly connected to ring atoms.
    sugar_region = set(sugar_ring_indices)
    for idx in sugar_ring_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # typically the 5'-CH2 group is a carbon; add if not already in the sugar ring.
            if nbr.GetIdx() not in sugar_region and nbr.GetSymbol() == 'C':
                sugar_region.add(nbr.GetIdx())
    
    # Step 2: Verify that a phosphate group is attached to the sugar region.
    # Look for any phosphorus which, via one oxygen neighbor, is connected to any atom in the sugar region.
    phosphate_attached = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # Phosphorus
            # iterate over neighbors of phosphorus
            for nbr in atom.GetNeighbors():
                # Looking for an oxygen connecting P to sugar:
                if nbr.GetSymbol() == 'O':
                    # Check if this oxygen is attached to any atom in the sugar region.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() in sugar_region:
                            phosphate_attached = True
                            break
                if phosphate_attached:
                    break
            if phosphate_attached:
                break
    if not phosphate_attached:
        return False, "No phosphate group attached to the sugar region via an oxygen bridge"
    
    # Step 3: Identify a nucleobase-like aromatic heterocycle.
    # Look for an aromatic ring (outside of our sugar region) that contains at least 2 nitrogen atoms.
    nucleobase_found = False
    for ring in ring_info.AtomRings():
        ring_set = set(ring)
        # Skip rings that are fully part of the sugar region.
        if ring_set.issubset(sugar_region):
            continue
        # Check aromaticity of all atoms in this ring.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        # Count nitrogen atoms in the ring.
        nitrogen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if nitrogen_count < 2:
            continue
        # Now check if at least one atom in the ring is directly attached to any atom in the sugar region.
        attached = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in sugar_region:
                    attached = True
                    break
            if attached:
                break
        if attached:
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "No nucleobase-like aromatic heterocycle (ring with ≥2 nitrogens attached to sugar) found"
    
    return True, "Molecule contains a nucleoside (sugar + nucleobase) with a phosphate group attached at a sugar hydroxy group"

# Example usage (for testing purposes; remove or modify if integrating elsewhere)
if __name__ == "__main__":
    # Testing with a few SMILES strings (from provided examples)
    test_smiles = [
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N",  # 2'-deoxy-5-methyl-5'-cytidylic acid (should be True)
        "[C@@H]1(N2C=C(C(=N)C=C2)C(O)=O)O[C@H](COP(O)(O)=O)[C@H]([C@H]1O)O",  # clitidine 5'-phosphate (should be True)
    ]
    for smi in test_smiles:
        result, reason = is_nucleotide(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")