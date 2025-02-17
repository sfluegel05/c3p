"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is defined as a carbohydrate derivative in which one or more of the
    oxygens or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR.
    
    The approach has two parts:
      1) It inspects ring systems of size 5 or 6 that have the expected characteristics
         for a sugar ring (mostly sp3 atoms with one heteroatom in the ring). We then check:
           A. If the intrinsic ring heteroatom is S,
           B. Or if any carbon in the ring bears an S substituent (unless that S is part of an
              extended (e.g. glucosinolate‐like) group).
      2) If no ring candidate is detected, the code performs an open‐chain search where we try to
         find a contiguous carbon chain (4–6 carbons long) that bears both hydroxyls and at least one
         thio substitution.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a thiosugar, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    sugar_candidate_found = False

    # Process each ring.
    for ring in ring_info.AtomRings():
        # Consider only rings of typical sugar sizes: 5 (furanose) or 6 (pyranose)
        if len(ring) not in (5, 6):
            continue

        # Identify the atoms that belong to the ring.
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        # For a typical sugar ring, we expect exactly one intrinsic heteroatom in the ring.
        ring_hetero_idxs = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
        if len(ring_hetero_idxs) != 1:
            continue

        # Check that all atoms in the ring are sp3 hybridized (sugar-like) 
        # and that several ring atoms carry external heteroatom substituents.
        valid_ring = True
        ext_het_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetHybridization() != rdchem.HybridizationType.SP3:
                valid_ring = False
                break
            # Count substituents (neighbors NOT in the ring) that are heteroatoms (O, N, S)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in ring and nbr.GetAtomicNum() in (7, 8, 16):
                    ext_het_count += 1
        if not valid_ring or ext_het_count < 2:
            continue
            
        # We have now a candidate sugar ring.
        sugar_candidate_found = True
        
        # CASE A: If the intrinsic heteroatom (in the ring) is sulfur, classify as thiosugar.
        ring_hetero = mol.GetAtomWithIdx(ring_hetero_idxs[0])
        if ring_hetero.GetSymbol() == "S":
            return True, "Sugar ring contains sulfur replacing the ring oxygen"
        
        # CASE B: Check if any carbon in the ring has an external substituent that is sulfur.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Only consider ring carbons.
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                # Skip atoms that are part of the ring.
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetSymbol() == "S":
                    # Check if the S appears to belong to an extended substituent (glucosinolate-like).
                    # We require that S be “simple” (e.g. not having extra heavy atom neighbors).
                    heavy_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() > 1 and a.GetIdx() != atom.GetIdx()]
                    if len(heavy_neighbors) <= 1:
                        return True, "Sugar ring has thio substitution at a hydroxyl-bearing carbon"
                    # Alternatively, check for double-bonded oxygen on the sulfur (extended functionality).
                    glucosinolate_like = False
                    for bond in nbr.GetBonds():
                        if bond.GetBondType() == rdchem.BondType.DOUBLE:
                            other = bond.GetOtherAtom(nbr)
                            if other.GetSymbol() == "O":
                                glucosinolate_like = True
                                break
                    if not glucosinolate_like:
                        return True, "Sugar ring has thio substitution at a hydroxyl-bearing carbon"
    # End ring processing.
    
    # If no sugar ring candidate was identified, try the open-chain search.
    # Here we look for a contiguous carbon chain (4–6 carbons) that bears several hydroxyl groups and one thio substitution.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        chain = []
        visited = set()
        # A simple depth-first search to get connected carbon chains (up to 6 atoms).
        def dfs(a, depth):
            if depth > 6:
                return
            chain.append(a.GetIdx())
            visited.add(a.GetIdx())
            for nbr in a.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                    dfs(nbr, depth + 1)
        dfs(atom, 1)
        # Check chain length typical of a monosaccharide.
        if len(chain) < 4 or len(chain) > 6:
            continue
        hydroxyl_count = 0
        thio_count = 0
        for idx in chain:
            a = mol.GetAtomWithIdx(idx)
            for nbr in a.GetNeighbors():
                if nbr.GetIdx() in chain:
                    continue
                if nbr.GetSymbol() == "O":
                    hydroxyl_count += 1
                if nbr.GetSymbol() == "S":
                    thio_count += 1
        # Require at least two hydroxyls (to indicate a carbohydrate) and one S.
        if hydroxyl_count >= 2 and thio_count >= 1:
            return True, "Open-chain carbohydrate derivative with thio substitution"
    
    # Final outcome: if a sugar ring was seen, but no thio substitution detected, or no sugar-like pattern was observed.
    if sugar_candidate_found:
        return False, "Sugar rings detected but no recognizable thio substitution found"
    else:
        return False, "No sugar-like pattern detected"

# Example usage (for testing, uncomment below)
#if __name__ == '__main__':
#    examples = {
#        "butylglucosinolic acid": "[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS(O)(=O)=O)/CCCC",
#        "1-(methylthio)ribulose 5-phosphate": "CSCC(=O)[C@H](O)[C@H](O)COP(O)(O)=O",
#        "6-thio-beta-D-galactose": "[C@@H]1([C@@H]([C@H]([C@H]([C@H](O1)CS)O)O)O)O"
#    }
#    for name, smi in examples.items():
#        result, reason = is_thiosugar(smi)
#        print(f"{name}: {result} ({reason})")