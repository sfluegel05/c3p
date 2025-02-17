"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: Acyl-CoA
Definition: An acyl-CoA is a thioester that results from the condensation of the thiol group 
of Coenzyme A with the carboxy group of a fatty acid. Its structure features an acyl thioester 
group directly attached to a pantetheine moiety (represented by the pattern SCCNC(=O)CCNC(=O)) 
which is connected to an adenosine nucleotide (detected via an adenine ring).
Improvement: Instead of only testing for the presence of the two fragments, we also check that 
the acyl thioester (CoA core) is topologically connected to the adenine motif via a short bond path.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    
    The function requires that the molecule contains:
      1. A thioester group directly connected to the CoA pantetheine fragment.
         This is captured by the SMARTS: "[CX3](=O)[SX2]CCNC(=O)CCNC(=O)".
      2. The adenine moiety commonly found in CoA should be present.
         This is detected by the SMARTS "n1cnc2c(N)ncnc12".
      3. In addition, at least one atom from the acyl-CoA core must be within a 
         short topological distance (<=15 bonds) of an atom in the adenine substructure.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an acyl-CoA, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Define the SMARTS for the acyl thioester core connected to the CoA pantetheine fragment.
    # It captures: acyl carbonyl C(=O) bonded to S, then “CCNC(=O)CCNC(=O)”
    core_smarts = "[CX3](=O)[SX2]CCNC(=O)CCNC(=O)"
    core_pattern = Chem.MolFromSmarts(core_smarts)
    if core_pattern is None:
        return False, "Error creating acyl-CoA core SMARTS pattern"
    
    core_matches = mol.GetSubstructMatches(core_pattern)
    if not core_matches:
        return False, "Acyl-CoA core not found (required thioester linked to pantetheine missing)"
    
    # 2. Define the SMARTS pattern for the adenine moiety.
    adenine_smarts = "n1cnc2c(N)ncnc12"
    adenine_pattern = Chem.MolFromSmarts(adenine_smarts)
    if adenine_pattern is None:
        return False, "Error creating adenine SMARTS pattern"
    
    adenine_matches = mol.GetSubstructMatches(adenine_pattern)
    if not adenine_matches:
        return False, "Adenine moiety of CoA not detected"
    
    # 3. Check connectivity between the acyl-CoA core and the adenine moiety.
    # In CoA these fragments are connected via the rigid pantetheine and phosphate linker.
    # For each pair of matched atoms (one from the core and one from adenine),
    # we compute the shortest topological distance (number of bonds).
    # If at least one pair is close (using a threshold of 15 atoms in the path), we consider
    # the core and adenine as part of one connected CoA structure.
    threshold = 15  # maximum number of atoms in the shortest path (empirical value)
    found_connection = False
    for core_match in core_matches:
        for adenine_match in adenine_matches:
            # For each pair of atoms (one from core_match, one from adenine_match),
            # compute the shortest path. We require at least one close connection.
            for core_atom_idx in core_match:
                for adenine_atom_idx in adenine_match:
                    path = Chem.GetShortestPath(mol, core_atom_idx, adenine_atom_idx)
                    if path and len(path) <= threshold:
                        found_connection = True
                        break
                if found_connection:
                    break
            if found_connection:
                break
        if found_connection:
            break
    
    if not found_connection:
        return False, "Acyl-CoA motifs found but not connected as expected"
    
    # Passed all tests: acyl thioester core exists, adenine is present, and they are connected.
    return True, "Contains acyl thioester linked to CoA pantetheine and adenine moiety"

# Example usage:
# status, reason = is_acyl_CoA("SMILES_STRING_HERE")
# print(status, reason)