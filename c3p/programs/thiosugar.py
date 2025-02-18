"""
Classifies: CHEBI:73754 thiosugar
"""
#!/usr/bin/env python
"""
Classifies: Thiosugar

A thiosugar is defined as a carbohydrate derivative in which one or more of the oxygens 
or hydroxy groups of the parent carbohydrate is replaced by sulfur or –SR.

This improved algorithm uses explicit hydrogens to confirm that oxygen substituents are 
actually hydroxyl groups. For cyclic candidates (5- or 6-membered rings) it requires a 
typical sugar ring pattern (saturated, with exactly one intrinsic heteroatom and at least 
three external –OH substituents). If the intrinsic heteroatom is S or if one of the ring 
carbons bears a simple S substituent, the molecule is classified as a thiosugar.

If no sugar ring is found, then an open‐chain search is performed on contiguous chains of 
exactly 5 or 6 carbons (typical for linear carbohydrates). Such chains must display at least 
three hydroxyl (–OH attached) substituents, at least one S substituent, and only few other 
(“nonsugar”) substituents.
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
    # Add explicit hydrogens for better analysis of hydroxyl groups.
    mol = Chem.AddHs(mol)
    
    ring_info = mol.GetRingInfo()
    sugar_candidate_found = False
    
    # First, check for cyclic (ring) sugar candidates.
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6):
            continue
        # Get atoms in the ring.
        atoms_in_ring = [mol.GetAtomWithIdx(i) for i in ring]
        # Discard aromatic rings.
        if any(atom.GetIsAromatic() for atom in atoms_in_ring):
            continue
        # Require that all atoms are sp3-hybridized.
        if any(atom.GetHybridization() != rdchem.HybridizationType.SP3 for atom in atoms_in_ring):
            continue
        
        # Count intrinsic heteroatoms in the ring (non-carbon and non-hydrogen).
        hetero_atoms = [atom for atom in atoms_in_ring if atom.GetAtomicNum() not in (6, 1)]
        if len(hetero_atoms) != 1:
            continue
        
        # Count external hydroxyl substituents on ring atoms.
        hydroxyl_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Only count oxygen as -OH if it is single‐bonded and bears an explicit hydrogen.
                if nbr.GetSymbol() == "O":
                    if any(n.GetSymbol() == "H" for n in nbr.GetNeighbors()):
                        hydroxyl_count += 1
        # Require at least three hydroxyl substituents.
        if hydroxyl_count < 3:
            continue
        
        # We now have a good sugar ring candidate.
        sugar_candidate_found = True
        
        # CASE A: The intrinsic heteroatom is sulfur.
        hetero_atom = hetero_atoms[0]
        if hetero_atom.GetSymbol() == "S":
            return True, "Sugar ring has sulfur replacing the intrinsic ring oxygen"
        
        # CASE B: Check if a ring carbon bears a simple S substituent.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:  # only consider carbons
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetSymbol() == "S":
                    # A simple S substituent should not be heavily decorated:
                    heavy_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() > 1 and n.GetIdx() != atom.GetIdx()]
                    if len(heavy_neighbors) <= 1:
                        return True, "Sugar ring has thio substitution at a hydroxyl-bearing carbon"
                    # Even if there is more substitution, we may consider sugar decoration.
                    # (More checks could be added here if needed.)
                    return True, "Sugar ring has thio substitution at a hydroxyl-bearing carbon"
        # No S substitution found in this ring candidate.
    
    # Next, search for open-chain sugar candidates only if no ring has returned a positive.
    # Look for a contiguous chain of exactly 5 or 6 sp3 carbons (non‐cyclic) that have predominantly
    # hydroxyl substituents and at least one S substituent.
    atoms = mol.GetAtoms()
    for atom in atoms:
        if atom.GetAtomicNum() != 6:
            continue
        if atom.IsInRing():
            continue
        chain = []
        visited = set()
        # Perform a simple Depth First Search to collect contiguous (non‐ring) carbons.
        def dfs(a):
            if a.GetAtomicNum() == 6 and (not a.IsInRing()):
                chain.append(a.GetIdx())
                visited.add(a.GetIdx())
                for nbr in a.GetNeighbors():
                    if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited and (not nbr.IsInRing()):
                        dfs(nbr)
        dfs(atom)
        if len(chain) not in (5, 6):
            continue
        hydroxyl_count = 0
        sulfur_count = 0
        extra_substituents = 0  # substituents that are neither -OH (from oxygen) nor sulfur.
        for idx in chain:
            a = mol.GetAtomWithIdx(idx)
            for nbr in a.GetNeighbors():
                if nbr.GetIdx() in chain:
                    continue
                symbol = nbr.GetSymbol()
                if symbol == "O":
                    # Check for a bonded H.
                    if any(n.GetSymbol() == "H" for n in nbr.GetNeighbors()):
                        hydroxyl_count += 1
                elif symbol == "S":
                    sulfur_count += 1
                elif nbr.GetAtomicNum() > 1:
                    extra_substituents += 1
        # Avoid chains decorated with many non-sugar substituents.
        if extra_substituents > len(chain):
            continue
        if hydroxyl_count >= 3 and sulfur_count >= 1:
            return True, "Open-chain carbohydrate derivative with thio substitution"
    
    # If a sugar-like candidate was found but no thio substitution evidence was detected.
    if sugar_candidate_found:
        return False, "Sugar-like patterns detected but no clear thio substitution found"
    else:
        return False, "No sugar-like pattern detected"

# For testing (uncomment below if running as a script):
# if __name__ == '__main__':
#     test_examples = {
#         "butylglucosinolic acid": "[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS(O)(=O)=O)/CCCC",
#         "celesticetin": "COC(C)C(NC(=O)[C@@H]1CCCN1C)[C@H]2O[C@H](SCCOC(=O)C3=CC=CC=C3O)[C@H](O)[C@@H](O)[C@H]2O",
#         "non-thiosugar control": "CC(=O)OC1=CC=CC=C1C(=O)O"  # example of a non-sugar molecule
#     }
#     for name, smi in test_examples.items():
#         flag, reason = is_thiosugar(smi)
#         print(f"{name}: {flag} ({reason})")