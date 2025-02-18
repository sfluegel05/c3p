"""
Classifies: CHEBI:73754 thiosugar
"""
#!/usr/bin/env python
"""
Classifies: Thiosugar

A thiosugar is defined as a carbohydrate derivative in which one or more of the oxygens
or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR.

The algorithm first searches for sugar rings (5- or 6-membered rings) that have
(sugar-like) external hydroxyl substituents. If such a ring is found,
either (A) its intrinsic heteroatom is S or (B) one of its ring carbons bears a simple S substituent.
If no ring candidate is found, it does an open-chain search for a contiguous chain of carbons (4-6)
that has predominantly oxygen substituents and at least one sulfur substituent.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    
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

    # Process each ring that is 5 or 6 members.
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6):
            continue
        
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count ring heteroatoms (non-carbon) and note the index of first one.
        ring_hetero_idxs = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6]
        # A typical sugar ring should have exactly one intrinsic heteroatom.
        if len(ring_hetero_idxs) != 1:
            continue

        # Check that all ring atoms are sp3 hybridized.
        valid_ring = True
        # Count hydroxyl (or oxygen) substituents on ring atoms.
        hydroxyl_subs = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetHybridization() != rdchem.HybridizationType.SP3:
                valid_ring = False
                break
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Count an O substituent only if it is attached as -OH (i.e. no double bonds)
                if nbr.GetSymbol() == "O":
                    # Check if O is connected to a hydrogen, assuming implicit H counts.
                    hydroxyl_subs += 1
        if not valid_ring or hydroxyl_subs < 2:
            continue  # not sugar-like enough

        # We have a sugar ring candidate
        sugar_candidate_found = True
        
        # CASE A: If the intrinsic heteroatom (in the ring) is sulfur, classify as thiosugar.
        ring_hetero = mol.GetAtomWithIdx(ring_hetero_idxs[0])
        if ring_hetero.GetSymbol() == "S":
            return True, "Sugar ring has sulfur replacing the intrinsic ring oxygen"

        # CASE B: Check if any carbon in the ring has a simple sulfur substituent.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:  # only checking carbons
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetSymbol() == "S":
                    # Check that this sulfur substituent is "simple": 
                    # It should have at most one heavy-atom (atomic num > 1) neighbor besides the ring carbon.
                    heavy_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() > 1 and a.GetIdx() != atom.GetIdx()]
                    if len(heavy_neighbors) <= 1:
                        return True, "Sugar ring has thio substitution at a hydroxyl-bearing carbon"
                    # If more heavy neighbors, check for possible carbonyl (e.g. glucosinolate-like) decoration.
                    glucosinolate_like = False
                    for bond in nbr.GetBonds():
                        if bond.GetBondType() == rdchem.BondType.DOUBLE:
                            other_atom = bond.GetOtherAtom(nbr)
                            if other_atom.GetSymbol() == "O":
                                glucosinolate_like = True
                                break
                    if not glucosinolate_like:
                        return True, "Sugar ring has thio substitution at a hydroxyl-bearing carbon"
    
    # End ring processing.

    # If no proper sugar ring candidate was detected, try an open-chain search.
    # The idea is to find a contiguous carbon chain of length 4–6 with mostly sugar-like substituents.
    atoms = mol.GetAtoms()
    for atom in atoms:
        if atom.GetAtomicNum() != 6:
            continue
        chain = []
        visited = set()
        # Simple DFS to get contiguous chain of carbons up to length 6.
        def dfs(a, depth):
            if depth > 6:
                return
            chain.append(a.GetIdx())
            visited.add(a.GetIdx())
            for nbr in a.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                    dfs(nbr, depth+1)
        dfs(atom, 1)
        if len(chain) < 4 or len(chain) > 6:
            continue

        # Check substituents on atoms in chain. We require at least 2 oxygen substituents (hydroxyls)
        # and at least one sulfur substituent.
        hydroxyl_count = 0
        sulfur_count = 0
        extra_substituents = 0  # count substituents that are neither O nor S
        for idx in chain:
            a = mol.GetAtomWithIdx(idx)
            for nbr in a.GetNeighbors():
                if nbr.GetIdx() in chain:
                    continue
                symbol = nbr.GetSymbol()
                if symbol == "O":
                    hydroxyl_count += 1
                elif symbol == "S":
                    sulfur_count += 1
                elif nbr.GetAtomicNum() > 1:
                    extra_substituents += 1
        # To be considered sugar-like, the chain should not have too many non-sugar substituents.
        if extra_substituents > len(chain):  
            continue
        if hydroxyl_count >= 2 and sulfur_count >= 1:
            return True, "Open-chain carbohydrate derivative with thio substitution"
    
    # If a sugar candidate was found but no thio substitution was detected, then report accordingly.
    if sugar_candidate_found:
        return False, "Sugar-like patterns detected but no clear thio substitution found"
    else:
        return False, "No sugar-like pattern detected"


# Example usage (for testing, uncomment the block below):
# if __name__ == '__main__':
#     test_examples = {
#         "butylglucosinolic acid": "[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS(O)(=O)=O)/CCCC",
#         "celesticetin": "COC(C)C(NC(=O)[C@@H]1CCCN1C)[C@H]2O[C@H](SCCOC(=O)C3=CC=CC=C3O)[C@H](O)[C@@H](O)[C@H]2O",
#         "non-thiosugar control": "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin - not a sugar
#     }
#     for name, smi in test_examples.items():
#         result, reason = is_thiosugar(smi)
#         print(f"{name}: {result} ({reason})")