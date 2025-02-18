"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: Nucleotide
A nucleotide is defined as a nucleoside phosphate resulting from the condensation 
of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.
This function tests for:
   - The presence of at least one phosphorus atom (signaling a phosphate group)
   - A furanose sugar ring (5-membered ring with exactly 1 oxygen and 4 carbons)
   - A nucleobase-like aromatic heterocycle (an aromatic ring system containing ≥2 nitrogen atoms)
   - A phosphate group attached to a sugar carbon via an oxygen (bridging oxygen of the hydroxy group)
"""
from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    
    A nucleotide should contain a nucleoside portion (a sugar ring attached to a nucleobase)
    with a phosphate group esterified at the 3' or 5' hydroxy position.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a nucleotide, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the presence of phosphorus (P, atomic number 15)
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphorus atom found (missing phosphate group)"
    
    # 2. Identify a furanose sugar ring: a 5-membered ring with exactly 1 oxygen and 4 carbons.
    sugar_ring_found = False
    sugar_ring_indices = set()
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 5:  # 5-membered ring
            num_ox = 0
            num_c = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetSymbol() == 'O':
                    num_ox += 1
                elif atom.GetSymbol() == 'C':
                    num_c += 1
            if num_ox == 1 and num_c == 4:
                sugar_ring_found = True
                sugar_ring_indices = set(ring)
                break
    if not sugar_ring_found:
        return False, "No furanose sugar ring (5-membered ring with 1 oxygen and 4 carbons) found"
    
    # 3. Identify a nucleobase-like aromatic heterocycle.
    # Instead of only accepting a single ring, we iterate over rings (other than the sugar ring)
    # and require that at least one aromatic ring system has at least 2 nitrogen atoms.
    nucleobase_found = False
    for ring in ring_info.AtomRings():
        ring_set = set(ring)
        # Skip the sugar ring if the ring is identical or completely within it.
        if ring_set == sugar_ring_indices or ring_set.issubset(sugar_ring_indices):
            continue
        # Check if all atoms in the ring are aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        # Count nitrogen atoms in the ring.
        n_nitrogen = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_nitrogen >= 2:
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "No nucleobase-like aromatic heterocycle (ring with ≥2 nitrogens) found"
    
    # 4. Verify the phosphate group is attached to the sugar by a bridging oxygen.
    # We allow for the case where the phosphate is not directly attached to an oxygen 
    # in the sugar ring but to an exocyclic oxygen (attached to a sugar carbon).
    phosphate_attached = False
    # For each phosphorus atom, examine its neighboring oxygens.
    for p in p_atoms:
        for neighbor in p.GetNeighbors():
            if neighbor.GetAtomicNum() != 8:  # We only care about oxygen neighbors
                continue
            # For this oxygen (bridging candidate), check its neighbors (other than the phosphorus)
            for nbr in neighbor.GetNeighbors():
                if nbr.GetIdx() == p.GetIdx():
                    continue
                # If the neighbor is an oxygen atom that belongs to the sugar ring, it's a direct bond.
                # Alternatively, if the neighbor is a carbon that is part of the sugar ring,
                # then this oxygen is the hydroxy oxygen from the sugar.
                if nbr.GetSymbol() == 'O' and nbr.GetIdx() in sugar_ring_indices:
                    phosphate_attached = True
                    break
                if nbr.GetSymbol() == 'C' and nbr.GetIdx() in sugar_ring_indices:
                    phosphate_attached = True
                    break
            if phosphate_attached:
                break
        if phosphate_attached:
            break

    if not phosphate_attached:
        return False, "Phosphate group is not attached to the sugar ring via a bridging oxygen (or its corresponding sugar carbon)"
    
    return True, "Molecule contains a nucleoside (furanose sugar + nucleobase) with a phosphate group attached"

# Example usage (you can remove these lines before deployment):
if __name__ == "__main__":
    # Example SMILES for 2'-deoxy-5-methyl-5'-cytidylic acid
    smiles_examples = [
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N",
        "[C@@H]1(N2C=C(C(=N)C=C2)C(O)=O)O[C@H](COP(O)(O)=O)[C@H]([C@H]1O)O"
    ]
    for smi in smiles_examples:
        res, reason = is_nucleotide(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")