"""
Classifies: CHEBI:28892 ganglioside
"""
"""
Classifies: Ganglioside – A molecule composed of a glycosphingolipid (ceramide and oligosaccharide)
with one or more sialic acids linked on the sugar chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    Gangliosides are complex molecules containing a ceramide (sphingolipid)
    linked to an oligosaccharide that bears one or more sialic acid residues.
    
    To improve our detection compared to previous attempts we:
      1. Add explicit hydrogens before computing substructures and ring information.
      2. Look for an amide bond typical of a ceramide (using the pattern "C(=O)N").
      3. Check for carboxylic acid groups in either protonated ("C(=O)O") or deprotonated ("C(=O)[O-]") forms,
         to catch the sialic acid markers.
      4. Detect sugar rings (representing the oligosaccharide) by:
           a. Matching a typical pyranose ring SMARTS pattern.
           b. Scanning all rings (of size 5–7 that are nonaromatic) and counting those with at least 2 oxygen atoms.
         At least two sugar rings are required.
      5. Check for a long chain of aliphatic (nonaromatic sp3) carbons as a proxy for the ceramide’s lipid chain.
         We require a longest chain of at least 16 carbons.
      6. Ensure the overall molecular weight is above a threshold (700 Da) to help weed out small molecules.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a ganglioside, False otherwise.
        str: Explanation or reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for better ring perception
    mol = Chem.AddHs(mol)
    
    # Heuristic 1: Check for an amide bond typical of a ceramide ("C(=O)N")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found—missing ceramide signature"
    
    # Heuristic 2: Check for carboxylic acid groups (both protonated and deprotonated forms)
    ca_pattern1 = Chem.MolFromSmarts("C(=O)O")
    ca_pattern2 = Chem.MolFromSmarts("C(=O)[O-]")
    if not (mol.HasSubstructMatch(ca_pattern1) or mol.HasSubstructMatch(ca_pattern2)):
        return False, "No carboxyl group found—sialic acid marker appears to be missing"
    
    # Heuristic 3: Identify sugar rings (oligosaccharide portion)
    # Strategy (a): Use a SMARTS for a typical pyranose ring.
    sugar_smarts = "O1C(CO)C(O)C(O)C1O"  # This is a loose pattern for a pyranose ring.
    sugar_mol = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches_smarts = mol.GetSubstructMatches(sugar_mol)
    
    # Strategy (b): Scan all rings for those of size 5-7, nonaromatic, with at least 2 oxygens.
    ring_info = mol.GetRingInfo()
    sugar_ring_count = 0
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6, 7):
            continue
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxy_count >= 2:
            sugar_ring_count += 1

    # Combine the two approaches for sugar detection
    total_sugar_count = max(len(sugar_matches_smarts), sugar_ring_count)
    if total_sugar_count < 2:
        return False, f"Only {total_sugar_count} sugar-like ring(s) detected; oligosaccharide portion appears insufficient"
    
    # Heuristic 4: Check for a long aliphatic (nonaromatic sp3 carbon) chain, as a proxy for the ceramide's fatty acid chain.
    longest_chain = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetIsAromatic() or atom.GetHybridization().name != "SP3":
            continue
        chain_length = 1
        visited = {atom.GetIdx()}
        stack = [(atom, chain_length)]
        while stack:
            current, length = stack.pop()
            for nbr in current.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (nbr.GetIdx() not in visited) and (nbr.GetHybridization().name == "SP3"):
                    visited.add(nbr.GetIdx())
                    new_len = length + 1
                    if new_len > longest_chain:
                        longest_chain = new_len
                    stack.append((nbr, new_len))
    if longest_chain < 16:
        return False, f"Longest aliphatic chain is too short ({longest_chain} carbons) for the ceramide part"
    
    # Heuristic 5: Check the overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical ganglioside"
    
    return True, "Contains the ceramide, oligosaccharide, and sialic acid features characteristic of a ganglioside"

# Example usage (uncomment below lines to test):
# result, reason = is_ganglioside("CCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)...")
# print(result, reason)