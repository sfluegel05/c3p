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

    Heuristics used in this implementation:
      1. Convert SMILES and add explicit hydrogens.
      2. Look for an amide bond ("C(=O)N") as a proxy for a ceramide portion.
      3. Look for at least one carboxyl group (C(=O)O or C(=O)[O-]) indicating a sialic acid (or similar).
      4. Search for sugar rings. Instead of the original 50% oxygen threshold,
         we now check nonaromatic rings of size 5–7 and call a ring “sugar‐like” if the fraction
         of oxygen atoms is at least 0.2. We require at least two such rings to flag an oligosaccharide.
      5. Assess the presence of a long aliphatic chain (nonaromatic sp3 carbon chain) of at least 16 carbons,
         to represent the lipid (fatty acid) portion of the ceramide.
      6. Require that the overall molecular weight is at least 700 Da.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is classified as a ganglioside, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to help with perception.
    mol = Chem.AddHs(mol)
    
    # Heuristic 1: Look for an amide bond as a probe for the ceramide portion.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found—missing ceramide signature"
    
    # Heuristic 2: Look for carboxyl groups as a proxy for sialic acids.
    ca_pattern1 = Chem.MolFromSmarts("C(=O)O")
    ca_pattern2 = Chem.MolFromSmarts("C(=O)[O-]")
    if not (mol.HasSubstructMatch(ca_pattern1) or mol.HasSubstructMatch(ca_pattern2)):
        return False, "No carboxyl group found—sialic acid marker appears to be missing"
    
    # Heuristic 3: Identify sugar rings (oligosaccharide part).
    ring_info = mol.GetRingInfo()
    sugar_ring_count = 0
    # For each ring, check size and count oxygen atoms. We now lower the fraction threshold to 0.2.
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6, 7):
            continue
        # Skip aromatic rings.
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if (oxy_count / len(ring)) >= 0.2:
            sugar_ring_count += 1
    if sugar_ring_count < 2:
        return False, f"Only {sugar_ring_count} sugar-like ring(s) detected; oligosaccharide portion appears insufficient"
    
    # Heuristic 4: Check for a long aliphatic chain representing the fatty acid tail.
    longest_chain = 0
    # Iterate over all carbon atoms that are sp3 and nonaromatic.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # not carbon
        if atom.GetIsAromatic() or atom.GetHybridization().name != "SP3":
            continue
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
    
    # Heuristic 5: Check overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical ganglioside"
    
    return True, "Contains ceramide, oligosaccharide and sialic acid features characteristic of a ganglioside"

# Example usage:
# result, reason = is_ganglioside("CCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)...")
# print(result, reason)