"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: Nucleotide
A nucleotide is defined as a nucleoside phosphate resulting from the condensation 
of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.
This improved function verifies:
   - Presence of a furanose sugar ring (5-membered ring with 1 oxygen and 4 carbons).
   - A phosphate ester attached to one of the sugar carbons via an oxygen.
   - A nucleobase-like aromatic heterocycle (an aromatic ring with ≥2 nitrogen atoms)
     that is bonded to the sugar.
"""
from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    
    A nucleotide should contain a nucleoside portion (a sugar ring attached to a nucleobase)
    with a phosphate group esterified at one of the sugar’s hydroxy positions.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a nucleotide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Identify a furanose sugar ring: a 5-membered ring with exactly 1 oxygen and 4 carbons.
    sugar_ring_found = False
    sugar_ring_indices = None
    ring_info = mol.GetRingInfo()
    # Iterate over all rings in the molecule.
    for ring in ring_info.AtomRings():
        if len(ring) == 5:  # furanose rings are 5-membered
            symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
            if symbols.count('O') == 1 and symbols.count('C') == 4:
                sugar_ring_found = True
                sugar_ring_indices = set(ring)
                break
    if not sugar_ring_found:
        return False, "No furanose sugar ring (5-membered ring with 1 oxygen and 4 carbons) found"
    
    # 2. Verify that a phosphate group is attached to the sugar via a bridging oxygen.
    # For each carbon atom in the sugar ring, check its neighbors for an oxygen that is bonded to a phosphorus.
    phosphate_attached = False
    for idx in sugar_ring_indices:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() != 'C':
            continue  # check only sugar carbons
        for nbr in atom.GetNeighbors():
            # Look for an oxygen substituent
            if nbr.GetSymbol() == 'O':
                # Check if this oxygen has a phosphorus neighbor.
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetAtomicNum() == 15:
                        phosphate_attached = True
                        break
                if phosphate_attached:
                    break
        if phosphate_attached:
            break
    if not phosphate_attached:
        return False, "No phosphate group attached to a sugar carbon via a bridging oxygen"
    
    # 3. Identify a nucleobase-like aromatic heterocycle.
    # Search for an aromatic ring (other than the sugar) that contains at least 2 nitrogen atoms and is connected to the sugar.
    nucleobase_found = False
    for ring in ring_info.AtomRings():
        ring_set = set(ring)
        # Skip the sugar ring.
        if sugar_ring_indices and ring_set & sugar_ring_indices == ring_set:
            continue
        # Ensure the ring is completely aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        # Count nitrogen atoms.
        n_nitrogens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if n_nitrogens >= 2:
            # Check if at least one atom in the aromatic ring is bonded to an atom in the sugar ring.
            attached = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in sugar_ring_indices:
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

# Example usage (for testing purposes; remove or modify if integrating elsewhere):
if __name__ == "__main__":
    # Example SMILES strings for nucleotides
    test_smiles = [
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N",  # 2'-deoxy-5-methyl-5'-cytidylic acid (expected True)
        "[C@@H]1(N2C=C(C(=N)C=C2)C(O)=O)O[C@H](COP(O)(O)=O)[C@H]([C@H]1O)O",  # clitidine 5'-phosphate (expected may fail nucleobase connectivity)
    ]
    for smi in test_smiles:
        result, reason = is_nucleotide(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")