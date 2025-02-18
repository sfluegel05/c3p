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

    Heuristics used:
      1. Add explicit hydrogens for better ring perception.
      2. Require an amide bond ("C(=O)N") as a proxy for the ceramide part.
      3. Require at least one carboxyl group (either "C(=O)O" or "C(=O)[O-]") which is common in sialic acids.
      4. Detect sugar rings by scanning all nonaromatic rings of size 5–7 and counting them
         as “sugar‐like” if the fraction of oxygen atoms in the ring is at least 0.5.
         (At least two such rings are required to confirm the oligosaccharide part.)
      5. Require a long aliphatic chain (nonaromatic sp3 carbon chain) of at least 16 carbons,
         a proxy for the lipid tail.
      6. Set an overall molecular weight threshold (at least 700 Da).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a ganglioside, False otherwise.
        str: Explanation or reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for better perception.
    mol = Chem.AddHs(mol)
    
    # Heuristic 1: Look for an amide bond ("C(=O)N") as a fingerprint for a ceramide.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found—missing ceramide signature"
    
    # Heuristic 2: Look for carboxylic acid group(s) (protonated or deprotonated) as proxy for sialic acid.
    ca_pattern1 = Chem.MolFromSmarts("C(=O)O")
    ca_pattern2 = Chem.MolFromSmarts("C(=O)[O-]")
    if not (mol.HasSubstructMatch(ca_pattern1) or mol.HasSubstructMatch(ca_pattern2)):
        return False, "No carboxyl group found—sialic acid marker appears to be missing"
    
    # Heuristic 3: Identify sugar rings (oligosaccharide portion).
    # Instead of a fixed SMARTS, we scan all rings that are nonaromatic with size 5-7.
    # We count a ring as sugar-like if the fraction of oxygen atoms is at least 0.5.
    ring_info = mol.GetRingInfo()
    sugar_ring_count = 0
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6, 7):
            continue
        # Skip rings that are aromatic.
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        # Count oxygen atoms in the ring.
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxy_count / len(ring) >= 0.5:
            sugar_ring_count += 1

    if sugar_ring_count < 2:
        return False, f"Only {sugar_ring_count} sugar-like ring(s) detected; oligosaccharide portion appears insufficient"
    
    # Heuristic 4: Check for a long aliphatic (nonaromatic sp3 carbon) chain to represent the ceramide fatty acid chain.
    longest_chain = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # Only carbon
        if atom.GetIsAromatic() or atom.GetHybridization().name != "SP3":
            continue  # Must be a nonaromatic sp3 carbon
        chain_length = 1
        visited = {atom.GetIdx()}
        stack = [(atom, chain_length)]
        while stack:
            current, length = stack.pop()
            for nbr in current.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and \
                   (nbr.GetIdx() not in visited) and (nbr.GetHybridization().name == "SP3"):
                    visited.add(nbr.GetIdx())
                    new_len = length + 1
                    if new_len > longest_chain:
                        longest_chain = new_len
                    stack.append((nbr, new_len))
    if longest_chain < 16:
        return False, f"Longest aliphatic chain is too short ({longest_chain} carbons) for the ceramide part"
    
    # Heuristic 5: Overall molecular weight check.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical ganglioside"
    
    return True, "Contains the ceramide, oligosaccharide and sialic acid features characteristic of a ganglioside"

# Example usage:
# result, reason = is_ganglioside("CCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)...")
# print(result, reason)