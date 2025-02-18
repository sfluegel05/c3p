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
    
    The following heuristic tests are applied:
      1. Check for an amide group (C(=O)N) typical of the ceramide part.
      2. Check for a carboxylic acid group (C(=O)[O;H]) that is often present in sialic acids.
      3. Identify sugar rings: we count rings having 5 or 6 atoms with at least 2 oxygen atoms.
         A glycosphingolipid should have an oligosaccharide (i.e. ≥2 sugar rings).
      4. Molecules that are gangliosides are usually large; we require a molecular weight above a threshold.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Heuristic 1: Check for amide bonds (expected in the ceramide part)
    # The SMARTS "C(=O)N" is used here as a minimal representation.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found—missing ceramide signature"
    
    # Heuristic 2: Look for carboxylic acid group (commonly present in sialic acids)
    ca_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    if not mol.HasSubstructMatch(ca_pattern):
        return False, "No carboxylic acid group found—sialic acid seems to be missing"
    
    # Heuristic 3: Identify sugar rings by scanning ring information.
    # We count 5- or 6-membered rings that contain at least 2 oxygen atoms.
    ring_info = mol.GetRingInfo()
    sugar_ring_count = 0
    for ring in ring_info.AtomRings():
        if len(ring) in (5, 6):
            # Count how many atoms in the ring are oxygen
            oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxy_count >= 2:
                sugar_ring_count += 1
    if sugar_ring_count < 2:
        return False, f"Only {sugar_ring_count} sugar-like ring(s) detected; oligosaccharide portion appears insufficient"
    
    # Heuristic 4: Check overall molecular weight.
    # Gangliosides are typically large molecules. (This cutoff is empirical and heuristic.)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, f"Molecular weight is too low ({mol_wt:.1f} Da) for a typical ganglioside"
    
    return True, "Contains ceramide, oligosaccharide, and sialic acid features characteristic of a ganglioside"

# Example usage:
# result, reason = is_ganglioside("CCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)...")
# print(result, reason)