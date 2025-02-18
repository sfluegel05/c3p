"""
Classifies: CHEBI:4194 D-hexose
"""
"""
Classifies: D-hexose
A D-hexose (for our purposes) is a hexose (exactly 6 carbons and 6 oxygens) that contains a sugar ring
(either a 5-membered or 6-membered ring with exactly one oxygen) and one ring-carbon that carries an
exocyclic –CH2OH substituent. The exocyclic CH2OH-bearing ring-carbon (the “C5” center) must have CIP
configuration 'R'.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_hexose(smiles: str):
    """
    Determines whether a molecule is a D-hexose.
    A D-hexose is defined here as a molecule that:
      - Has exactly 6 carbon atoms and 6 oxygen atoms.
      - Contains a sugar ring (a 5- or 6-membered ring containing exactly 1 oxygen atom).
      - Has exactly one exocyclic CH2OH group attached to a ring carbon. This ring carbon (the C5 center)
        must have the CIP configuration 'R'.
        
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        (bool, str): Tuple with classification result and a message explaining the decision.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Assign stereochemistry so that CIP codes get computed.
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    # Check for exactly 6 carbons and 6 oxygens.
    atoms = list(mol.GetAtoms())
    carbons = [atom for atom in atoms if atom.GetAtomicNum() == 6]
    oxygens = [atom for atom in atoms if atom.GetAtomicNum() == 8]
    if len(carbons) != 6:
        return False, "Molecule does not have exactly 6 carbon atoms (not a hexose)"
    if len(oxygens) != 6:
        return False, "Molecule does not have exactly 6 oxygen atoms (likely a modified hexose)"
    
    # Look for an exocyclic CH2OH group.
    # In a proper CH2OH group the carbon should be:
    #  - sp3 carbon (implicit in rdkit) with exactly 2 hydrogens,
    #  - bonded to exactly two heavy atoms (one oxygen and one ring carbon).
    exo_candidates = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # We require that this carbon is NOT in a ring.
        if atom.IsInRing():
            continue
        # Check the number of attached hydrogens.
        if atom.GetTotalNumHs() != 2:
            continue
        # Count heavy neighbors (atoms with atomic number > 1)
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 2:
            continue
        # We expect one neighbor to be oxygen (for the –OH) and the other to be a carbon.
        nbr_atomic_nums = sorted([nbr.GetAtomicNum() for nbr in heavy_neighbors])
        if nbr_atomic_nums != [6, 8]:
            continue
        # Identify the neighbor that is a carbon.
        ring_carbon = None
        for nbr in heavy_neighbors:
            if nbr.GetAtomicNum() == 6:
                # For a proper exocyclic attachment, the carbon neighbor should be in a ring.
                if nbr.IsInRing():
                    ring_carbon = nbr
        if ring_carbon is None:
            continue
        # This atom qualifies as a CH2OH substituent attached to a ring carbon.
        exo_candidates.append((atom, ring_carbon))
    
    if len(exo_candidates) == 0:
        return False, "No appropriate exocyclic CH2OH substituent found"
    if len(exo_candidates) > 1:
        return False, "Multiple exocyclic CH2OH substituents found (ambiguous candidate for C5)"
    
    # We have exactly one exocyclic CH2OH candidate.
    exo_atom, candidate_C5 = exo_candidates[0]
    
    # Check that candidate_C5 is in a sugar ring.
    # A typical sugar ring is either 5-membered (furanose) or 6-membered (pyranose)
    # and must contain exactly one oxygen.
    ring_info = mol.GetRingInfo()
    sugar_ring = None
    for ring in ring_info.AtomRings():
        if len(ring) not in (5, 6):
            continue
        if candidate_C5.GetIdx() not in ring:
            continue
        # Count oxygen atoms in the ring.
        oxygen_count = 0
        for idx in ring:
            if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8:
                oxygen_count += 1
        if oxygen_count != 1:
            continue
        # Accept the first ring that meets the criteria.
        sugar_ring = ring
        break
    
    if sugar_ring is None:
        return False, "Candidate C5 is not part of an appropriate sugar ring (5- or 6-membered with one oxygen)"
    
    # Now check the stereochemistry of the candidate C5 atom.
    # The CIP code should be available (R expected for D-hexose).
    if not candidate_C5.HasProp('_CIPCode'):
        return False, "Stereochemistry was not assigned for the candidate C5 atom"
    cip = candidate_C5.GetProp('_CIPCode')
    if cip != 'R':
        return False, f"C5 has CIP configuration {cip} (expected R for D-hexose)"
    
    return True, "Molecule is a D-hexose: 6 carbons with a sugar ring and the C5 center in R configuration"


# Example usage:
if __name__ == "__main__":
    # Test with a couple of SMILES strings for D-hexoses.
    test_smiles_list = [
        "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@H]1O",  # alpha-D-galactose
        "O1[C@@H]([C@@H](O)[C@H](O)[C@@H]1O)[C@H](O)CO"   # beta-D-idopyranose
    ]
    for smi in test_smiles_list:
        res, reason = is_D_hexose(smi)
        print(f"SMILES: {smi}\nResult: {res}\nReason: {reason}\n{'-'*40}")