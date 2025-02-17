"""
Classifies: CHEBI:36976 nucleotide
"""
"""
Classifies: Nucleotide
A nucleotide is defined as a nucleoside phosphate resulting from the condensation 
of the 3 or 5 hydroxy group of a nucleoside with phosphoric acid.

This function checks that:
  - A furanose sugar ring (5-membered ring with exactly 1 oxygen and 4 carbons) is present.
  - The sugar region is expanded with its immediate carbon neighbors (e.g. the exocyclic 5'-CH2 group).
  - A glycosidic bond exists between a carbon atom of the sugar (typically the anomeric carbon) 
    and a nucleobase-like fragment (an aromatic heterocycle carrying at least one nitrogen).
  - A phosphate group is attached to the sugar through an oxygen (i.e. a P atom with an O that connects to the sugar region).
  
If these criteria are met, the molecule is classified as a nucleotide.
"""
from rdkit import Chem

def is_nucleotide(smiles: str):
    """
    Determines if a molecule is a nucleotide based on its SMILES string.
    
    A nucleotide is expected to have:
      • A nucleoside portion: a furanose sugar (and its immediate exocyclic carbons) 
        connected via a glycosidic bond to a nucleobase (an aromatic heterocycle with at least one nitrogen).
      • A phosphate group attached to the sugar via an oxygen bridge.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a nucleotide, False otherwise.
        str: Reason for the classification result.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    sugar_ring = None
    # Step 1: Identify a furanose ring – a 5-membered ring with exactly 1 oxygen and 4 carbons.
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            # Count atoms in the ring by element.
            symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
            if symbols.count("O") == 1 and symbols.count("C") == 4:
                sugar_ring = set(ring)
                break
    if sugar_ring is None:
        return False, "No furanose sugar ring (5-membered ring with 1 oxygen and 4 carbons) found"
    
    # Step 2: Expand sugar region: include immediate carbon neighbors of the sugar ring.
    sugar_region = set(sugar_ring)
    for idx in list(sugar_ring):
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # Consider adding carbons that lie immediately outside the ring.
            if nbr.GetSymbol() == "C":
                sugar_region.add(nbr.GetIdx())
    
    # Step 3: Identify the glycosidic bond.
    # Look for a bond from a carbon atom in the sugar ring to a nitrogen that belongs to an aromatic ring (nucleobase).
    glycosidic_found = False
    for idx in sugar_ring:
        sugar_atom = mol.GetAtomWithIdx(idx)
        # we expect the connecting atom (anomeric carbon) to be carbon
        if sugar_atom.GetSymbol() != "C":
            continue
        for nbr in sugar_atom.GetNeighbors():
            if nbr.GetIdx() in sugar_region:
                # Skip atoms that are also part of (or immediately attached to) the sugar
                continue
            # We favor a connection to a nitrogen.
            if nbr.GetAtomicNum() == 7 and nbr.GetIsAromatic():
                # Ensure that the neighbor is part of an aromatic ring that is not completely in the sugar region.
                for ring in ring_info.AtomRings():
                    if nbr.GetIdx() in ring and not set(ring).issubset(sugar_region):
                        glycosidic_found = True
                        break
            if glycosidic_found:
                break
        if glycosidic_found:
            break
    if not glycosidic_found:
        return False, ("No glycosidic bond found connecting a sugar carbon to an aromatic heterocycle "
                       "with nitrogen (nucleobase-like component)")
    
    # Step 4: Check for a phosphate group attached to the sugar.
    # We consider a phosphate as a phosphorus atom (atomic number 15) that is connected 
    # via an oxygen to any atom in the sugar region.
    phosphate_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # phosphorus
            for nbr in atom.GetNeighbors():
                if nbr.GetSymbol() == "O":
                    # Check if this oxygen is connected to any atom in the sugar region.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() in sugar_region:
                            phosphate_found = True
                            break
                if phosphate_found:
                    break
        if phosphate_found:
            break
    if not phosphate_found:
        return False, "No phosphate group attached to the sugar region via an oxygen bridge"
    
    return True, ("Molecule contains a furanose sugar (with expanded sugar region), a glycosidic bond to a nucleobase-like "
                  "aromatic heterocycle, and a phosphate group attached via oxygen – consistent with a nucleotide")
    
# Example usage (for testing purposes)
if __name__ == "__main__":
    # List a few test examples; note that in practice the full test list is extensive.
    test_examples = [
        "Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)nc1N",  # 2'-deoxy-5-methyl-5'-cytidylic acid (should be True)
        "[C@@H]1(N2C=C(C(=N)C=C2)C(O)=O)O[C@H](COP(O)(O)=O)[C@H]([C@H]1O)O",  # clitidine 5'-phosphate (should be True)
    ]
    for smi in test_examples:
        res, reason = is_nucleotide(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n")