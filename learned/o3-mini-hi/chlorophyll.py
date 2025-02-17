"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: Chlorophyll – a family of magnesium porphyrins defined by the presence 
of a fifth ring beyond the four pyrrole-like rings and usually a long phytol chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophyll is defined (approximately) as a magnesium porphyrin that meets the 
    following criteria:
      1. Contains at least one magnesium atom.
      2. The magnesium is coordinated to four nitrogen atoms.
      3. Contains at least five rings overall, and one of those rings is a five-membered ring
         (interpreted as the extra isocyclic ring in chlorophylls).
      4. Contains a long aliphatic chain (phytol chain) – here we search for an unbroken chain 
         of at least 10 carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as chlorophyll, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of magnesium (Mg, atomic number 12)
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if not mg_atoms:
        return False, "No magnesium atom found"
    
    # Verify that at least one Mg atom is coordinated to four nitrogen atoms.
    mg_with_four_n = False
    for mg in mg_atoms:
        n_neighbors = sum(1 for nbr in mg.GetNeighbors() if nbr.GetAtomicNum() == 7)
        if n_neighbors >= 4:
            mg_with_four_n = True
            break
    if not mg_with_four_n:
        return False, "Magnesium is not coordinated to four nitrogen atoms (porphyrin core missing)"
    
    # Check for the presence of a long phytol chain.
    # Here we define a simple aliphatic chain of 10 carbons as a proxy.
    phytol_pattern = Chem.MolFromSmarts("CCCCCCCCCC")
    if not mol.HasSubstructMatch(phytol_pattern):
        return False, "Long phytol chain (10 consecutive carbons) not found"
    
    # Analyze the ring structure. Chlorophyll typically has a porphyrin macrocycle (4 fused rings)
    # plus an additional (isocyclic) ring – thus at least 5 rings in the core.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 5:
        return False, f"Insufficient number of rings found ({num_rings}); expected at least 5"
    
    # Look for a five-membered ring (the extra ring typical in chlorophyll) among the rings.
    atom_rings = ring_info.AtomRings()
    has_five_membered = any(len(ring) == 5 for ring in atom_rings)
    if not has_five_membered:
        return False, "No five-membered ring (extra isocyclic ring) found"
    
    return True, "Contains magnesium porphyrin core (Mg coordinated by 4 N), extra five-membered ring and phytol chain"
    
# Example usage (you may comment out these lines when using the function in another module):
if __name__ == "__main__":
    test_smiles = "CCC1=C(C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31"  # Example: chlorophyll a
    result, reason = is_chlorophyll(test_smiles)
    print(result, reason)