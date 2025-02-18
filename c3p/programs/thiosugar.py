"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is defined as a carbohydrate derivative in which one or more
    of the oxygens or hydroxy groups of the parent carbohydrate are replaced by sulfur or -SR.
    
    The algorithm works in two parts:
      1) It inspects ring systems of size 5 or 6 that have the expected characteristics
         for a sugar ring (most atoms are sp3 carbon with one heteroatom, normally O).
         If the one heteroatom is S instead of O, or if any ring carbon that would normally
         bear a hydroxyl substituent instead carries a sulfur substituent (and that S is not
         part of an extended substituent such as in glucosinolates), then the molecule is
         classified as a thiosugar.
      2) If no sugar ring candidate is detected the code performs a crude open‐chain search:
         if a contiguous chain of 4–6 carbons shows evidence of thio substitution replacing
         an expected hydroxyl group, the molecule is flagged.
    
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
    
    ri = mol.GetRingInfo()
    sugar_found = False

    # Process each ring in the molecule.
    for ring in ri.AtomRings():
        # Consider only rings of typical sugar sizes: 5 (furanose) or 6 (pyranose)
        if len(ring) not in (5, 6):
            continue

        # Identify non-carbon atoms in the ring.
        non_c = [mol.GetAtomWithIdx(i) for i in ring if mol.GetAtomWithIdx(i).GetAtomicNum() != 6]
        # For a typical sugar ring, we expect exactly one heteroatom.
        if len(non_c) != 1:
            continue

        # Further check that the atoms in the ring have sugar-like properties: most should be sp3 hybridized,
        # and several should have external heteroatom substituents (e.g. -OH).
        sugar_candidates = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                break
            # Count substituents (neighbors outside the ring) that are heteroatoms (N, O, or S).
            ext_hets = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring and nbr.GetAtomicNum() in (7, 8, 16)]
            sugar_candidates.append(len(ext_hets))
        else:
            # If we did not break, require that a minimum number of ring atoms bear heteroatom substituents.
            if sum(sugar_candidates) < 2:
                continue
            sugar_found = True  # We have a candidate sugar ring.
            
            # CASE A: Check if the heteroatom on the ring is sulfur.
            hetero = non_c[0]
            if hetero.GetSymbol() == 'S':
                return True, "Sugar ring contains sulfur in place of oxygen"
            
            # CASE B: Check if any carbon in the ring has a sulfur substituent instead of a hydroxyl group.
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() != 6:
                    continue
                # Look at neighbors outside the ring.
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() in ring:
                        continue
                    if nbr.GetSymbol() == 'S':
                        # To avoid misclassification with glucosinolate-like groups (where sulfur is part of an O=S linkage),
                        # check if the sulfur has extra connectivities indicating extended functionality.
                        s_neighbors = [a for a in nbr.GetNeighbors() if a.GetIdx() != atom.GetIdx()]
                        glucosinolate_like = False
                        if s_neighbors:
                            for sn in s_neighbors:
                                # Examine bonds for a double bond to oxygen.
                                for bond in nbr.GetBonds():
                                    if (bond.GetBeginAtom() == nbr and bond.GetEndAtom().GetSymbol() == 'O' and 
                                        bond.GetBondTypeAsDouble() == 2.0):
                                        glucosinolate_like = True
                                        break
                                    if (bond.GetEndAtom() == nbr and bond.GetBeginAtom().GetSymbol() == 'O' and 
                                        bond.GetBondTypeAsDouble() == 2.0):
                                        glucosinolate_like = True
                                        break
                                if glucosinolate_like:
                                    break
                        if glucosinolate_like:
                            continue
                        return True, "Sugar ring has thio substitution at a hydroxyl-bearing carbon"
    # End ring processing

    # Open-chain search if no sugar ring candidate was detected.
    if not sugar_found:
        atoms = mol.GetAtoms()
        for atom in atoms:
            if atom.GetAtomicNum() != 6:
                continue
            chain = []
            visited = set()
            # Simple DFS to collect connected carbon chains (up to 6 atoms).
            def dfs(a, depth):
                if depth > 6:
                    return
                chain.append(a.GetIdx())
                visited.add(a.GetIdx())
                for nbr in a.GetNeighbors():
                    if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                        dfs(nbr, depth + 1)
                return
            dfs(atom, 1)
            if not (4 <= len(chain) <= 6):
                continue
            n_O = 0
            n_S = 0
            for idx in chain:
                a = mol.GetAtomWithIdx(idx)
                for nbr in a.GetNeighbors():
                    if nbr.GetIdx() in chain:
                        continue
                    if nbr.GetSymbol() == 'O':
                        n_O += 1
                    if nbr.GetSymbol() == 'S':
                        n_S += 1
            if n_S >= 1 and n_O >= 1:
                return True, "Open-chain carbohydrate derivative with thio substitution"
    
    # Final determination based on what was found.
    if sugar_found:
        return False, "Sugar rings detected but no thio substitution found"
    else:
        return False, "No sugar (ring or open-chain) pattern detected"
    
# (End of code)
    
# Example usage (uncomment to test):
#if __name__ == '__main__':
#    test_examples = {
#        "butylglucosinolic acid": "[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS(O)(=O)=O)/CCCC",
#        "1-(methylthio)ribulose 5-phosphate": "CSCC(=O)[C@H](O)[C@H](O)COP(O)(O)=O",
#        "6-thio-beta-D-galactose": "[C@@H]1([C@@H]([C@H]([C@H]([C@H](O1)CS)O)O)O)O"
#    }
#    for name, smi in test_examples.items():
#        result, reason = is_thiosugar(smi)
#        print(f"{name}: {result} ({reason})")