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
    Gangliosides are complex lipids characterized by:
      - A ceramide portion (represented here by the presence of at least one amide bond).
      - An oligosaccharide chain containing one or more sialic acid residues (indicated by at least one carboxyl group and sugar rings).
      - A long aliphatic chain (many sp3 carbon atoms) to represent the fatty acid part.
      - A relatively high molecular weight (here, at least 700 Da).
      
    Improvements from the previous approach:
      - Instead of requiring two sugar rings with a high oxygen fraction, we count any ring of size 5–7 (nonaromatic) 
        that contains at least one oxygen as a plausible sugar ring. We require at least one such ring.
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a ganglioside, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to improve substructure perception.
    mol = Chem.AddHs(mol)
    
    # Heuristic 1: Look for at least one amide bond as a proxy for the ceramide portion.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found – missing ceramide signature"
    
    # Heuristic 2: Look for carboxyl groups indicating sialic acid(s).
    # Check both deprotonated and neutral carboxyl forms.
    carboxyl1 = Chem.MolFromSmarts("C(=O)O")
    carboxyl2 = Chem.MolFromSmarts("C(=O)[O-]")
    if not (mol.HasSubstructMatch(carboxyl1) or mol.HasSubstructMatch(carboxyl2)):
        return False, "No carboxyl group found – sialic acid marker appears to be missing"
    
    # Heuristic 3: Identify sugar rings within the oligosaccharide portion.
    # We now count any nonaromatic ring of size 5, 6, or 7 that contains at least one oxygen atom.
    ring_info = mol.GetRingInfo()
    sugar_ring_count = 0
    for ring in ring_info.AtomRings():
        # Consider only rings of size 5-7.
        if len(ring) not in (5, 6, 7):
            continue
        # Skip the ring if any atom in it is aromatic.
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        # Count this ring as "sugar-like" if at least one atom is oxygen.
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring):
            sugar_ring_count += 1
    if sugar_ring_count < 1:
        return False, f"Only {sugar_ring_count} sugar-like ring(s) detected; oligosaccharide portion appears insufficient"
    
    # Heuristic 4: Check for a long aliphatic chain representing the fatty acid tail.
    longest_chain = 0
    # Iterate over all carbon atoms that are sp3 and nonaromatic.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # not a carbon atom
        if atom.GetIsAromatic() or atom.GetHybridization().name != "SP3":
            continue
        chain_length = 1
        visited = {atom.GetIdx()}
        stack = [(atom, chain_length)]
        while stack:
            current, length = stack.pop()
            for nbr in current.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and \
                   (nbr.GetHybridization().name == "SP3") and (nbr.GetIdx() not in visited):
                    visited.add(nbr.GetIdx())
                    new_length = length + 1
                    if new_length > longest_chain:
                        longest_chain = new_length
                    stack.append((nbr, new_length))
    if longest_chain < 16:
        return False, f"Longest aliphatic chain is too short ({longest_chain} carbons) for the ceramide part"
    
    # Heuristic 5: Check overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical ganglioside"
    
    return True, "Contains ceramide (amide bond), oligosaccharide (sugar-like ring and carboxyl group) and a long fatty acid chain characteristic of a ganglioside"

# Example usage:
# result, reason = is_ganglioside("CCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)...")
# print(result, reason)