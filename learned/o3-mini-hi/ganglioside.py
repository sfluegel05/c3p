"""
Classifies: CHEBI:28892 ganglioside
"""
We can see from the prior attempt that our simple checks (a substructure for an amide bond, one for a carboxylic acid, plus “counting” sugar‐rings by scanning the standard ring information) are giving us two problems. On one hand, a molecule like eliglustat tartrate (which is not a ganglioside) is matching the “has ceramide and carboxyl” tests even though it lacks a bona fide glycosphingolipid portion. On the other hand, many actual gangliosides (which by definition have one or more “sugar” rings in their oligosaccharide part) are not being flagged because our “sugar‐ring” detection (based on simply counting rings with ≥2 oxygen atoms inside 5–6 members) does not work when the sugars are “decorated” by multiple substituents or have stereochemistry that prevents the simple match from working. 

The improvements we can try are as follows:
 
1. When parsing the molecule, add explicit hydrogens so that our ring perception has the full valence information.
2. Instead of one single SMARTS for carboxylic acid we can look for either protonated or deprotonated forms (using two patterns) so that we are not fooled by ionization state.
3. For sugars (i.e. the oligosaccharide part) we use a two‐pronged strategy: (a) use a SMARTS for a “typical” pyranose ring (even if this is not perfect) and also (b) scan the molecule’s ring info (after adding hydrogens) for rings that are nonaromatic, of size 5–7, and that have at least two oxygen atoms. (A match by either procedure will count as a sugar ring.)
4. In addition, we add a check on the ceramide portion. Gangliosides always come from glycosphingolipids so there should be (i) a ceramide “signature” (an amide group joined to a long alkyl chain and a sphingoid base) and (ii) a long lipid chain present. We approximate this by requiring that the atom count of the longest alkane chain in the molecule is at least 16.
5. We also keep a molecular‐weight cutoff (which can be tweaked) to help weed out smaller molecules.

Below is a single Python code block that implements these ideas. (Note that one can improve on each heuristic; the approach below is a first approximation.) 

Below is the revised code:

------------------------------------------------

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
    
    To improve on our previous attempt we now:
      1. Add explicit hydrogens before computing substructures and ring information.
      2. Look for amide groups as part of the ceramide.
      3. Check for carboxylic acid groups (both protonated and deprotonated) which are common in sialic acids.
      4. Detect sugar rings using two strategies:
           a. A SMARTS based on a typical (non‐aromatic) pyranose ring.
           b. A general scan over rings (of size 5–7 and nonaromatic) counting rings with at least 2 oxygens.
         We require that at least two sugar units are present.
      5. Ensure there is a very long aliphatic chain (a proxy for the sphingoid/ceramide lipid chain).
      6. Check that the overall molecular weight is above an empirical threshold.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a ganglioside, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that ring perception is more reliable.
    mol = Chem.AddHs(mol)
    
    # Heuristic 1: Check for amide bond typical of ceramide.
    # Look for the substructure "C(=O)N"
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found—missing ceramide signature"
    
    # Heuristic 2: Check for carboxylic acid groups.
    # Try both protonated ("C(=O)O") and deprotonated ("C(=O)[O-]") forms.
    ca_pattern1 = Chem.MolFromSmarts("C(=O)O")
    ca_pattern2 = Chem.MolFromSmarts("C(=O)[O-]")
    if not (mol.HasSubstructMatch(ca_pattern1) or mol.HasSubstructMatch(ca_pattern2)):
        return False, "No carboxyl group found—sialic acid marker is missing"
    
    # Heuristic 3: Identify sugar rings (oligosaccharide portion).
    # Strategy (a): Use a SMARTS for a typical pyranose ring.
    sugar_smarts = "O1C(CO)C(O)C(O)C1O"  # a loose pattern for a pyranose ring (stereochemistry ignored)
    sugar_mol = Chem.MolFromSmarts(sugar_smarts)
    sugar_matches_smarts = mol.GetSubstructMatches(sugar_mol)
    
    # Strategy (b): Scan through each ring in the molecule, requiring size 5-7, nonaromatic, with >=2 oxygens.
    ring_info = mol.GetRingInfo()
    sugar_ring_count = 0
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6, 7):
            continue
        # check that the ring is nonaromatic
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
        if oxy_count >= 2:
            sugar_ring_count += 1

    # Combine sugar ring detections (count unique rings from either strategy)
    total_sugar_count = max(len(sugar_matches_smarts), sugar_ring_count)
    if total_sugar_count < 2:
        return False, f"Only {total_sugar_count} sugar-like ring(s) detected; oligosaccharide portion appears insufficient"

    # Heuristic 4: Check for a long aliphatic (lipid) chain as a proxy for the ceramide's fatty acid chain.
    # We approximate this by computing the length of the longest chain consisting only of aliphatic carbons.
    longest_chain = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Consider only non-aromatic sp3 carbons.
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
        return False, f"Longest aliphatic chain is too short ({longest_chain} carbons) for a ceramide part"
    
    # Heuristic 5: Check overall molecular weight.
    # Gangliosides are typically large molecules.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, f"Molecular weight is too low ({mol_wt:.1f} Da) for a typical ganglioside"
    
    return True, "Contains ceramide, oligosaccharide, and sialic acid features characteristic of a ganglioside"

# Example usage (uncomment to test):
# result, reason = is_ganglioside("CCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)...")
# print(result, reason)
  
"""
Notes:
- This heuristic approach is still approximate and may misclassify molecules.
- The sugar SMARTS and ring criteria were relaxed (and explicit hydrogens added)
  to try to capture more of the diverse glycosidic motifs.
- Likewise, the longest-chain algorithm is a rough proxy for a ceramide’s lipid chain.
- Further refinements (or using a specialized glycosphingolipid classifier) may be needed.
"""
