"""
Classifies: CHEBI:35681 secondary alcohol
"""
"""
Classifies: Secondary Alcohol
Definition:
  A secondary alcohol is defined here as a molecule in which every non–acidic –OH group 
  (after adding explicit hydrogens) is attached to a saturated carbon (sp³) that is 
  bound to exactly two carbons and one hydrogen. 
Note:
  (1) –OH groups that are part of a carboxylic acid (or its deprotonated form) are ignored.
  (2) This is a heuristic approach. Some complex molecules (e.g. steroids) may be mis‐classified.
"""
from rdkit import Chem

def is_secondary_alcohol(smiles: str):
    """
    Determines if a molecule qualifies as a secondary alcohol. In our implementation:
      - We first parse the SMILES and add explicit hydrogens.
      - We then ignore –OH groups that belong to carboxylic acids (or carboxylate forms).
      - For every remaining (candidate) –OH, we require that the oxygen is attached (via a single bond)
        to a saturated carbon atom (sp³) and that (aside from the –OH oxygen) that carbon is bonded to exactly 
        two carbons and one hydrogen.
      - Only if at least one candidate –OH is found and every candidate meets the secondary criteria do we classify
        the molecule as a secondary alcohol.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a secondary alcohol, False otherwise.
        str: Explanation/reason for the decision.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens.
    mol = Chem.AddHs(mol)

    # First, detect the –OH groups that are part of carboxylic acids.
    # Carboxylic acid OH groups have a SMARTS pattern like C(=O)[OX1H] (or with a negative charge).
    acid_oh_smarts = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    acid_oh_idxs = set()
    for match in mol.GetSubstructMatches(acid_oh_smarts):
        # In the pattern, the oxygen is the third atom.
        if len(match) >= 3:
            acid_oh_idxs.add(match[2])
    
    total_candidate = 0
    secondary_candidate = 0
    
    # Iterate over all atoms looking for oxygen atoms that have a hydrogen (i.e. are –OH) 
    # and are not part of a carboxylic acid (from our SMARTS).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        idx = atom.GetIdx()
        # Skip if this oxygen is already identified as part of a carboxylic acid.
        if idx in acid_oh_idxs:
            continue
        
        # For a non–acidic –OH, we expect at least one hydrogen attached.
        # (Sometimes an oxygen could be bridging; here we require exactly one hydrogen partner)
        neighbors = atom.GetNeighbors()
        # Get the indices of hydrogen and carbon neighbors.
        h_neighbors = [nb for nb in neighbors if nb.GetAtomicNum() == 1]
        c_neighbors = [nb for nb in neighbors if nb.GetAtomicNum() == 6]
        
        if len(h_neighbors) < 1 or len(c_neighbors) != 1:
            # If there isn’t exactly one carbon neighbor (or at least one hydrogen), skip.
            continue
        
        # We have found a candidate –OH group.
        total_candidate += 1
        # Let the single carbon neighbor be the one bearing the OH.
        carbon = c_neighbors[0]
        # Check that the carbon is sp³ (i.e. saturated).
        if carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue  # not on a saturated carbon
        
        # Now, examine the substituents on that carbon (excluding the –OH itself).
        count_carbons = 0
        count_hydrogens = 0
        for nb in carbon.GetNeighbors():
            if nb.GetIdx() == atom.GetIdx():
                continue  # skip the –OH oxygen
            if nb.GetAtomicNum() == 6:
                count_carbons += 1
            elif nb.GetAtomicNum() == 1:
                count_hydrogens += 1
        # For a secondary alcohol carbon, we want exactly two carbon substituents and one hydrogen.
        if count_carbons == 2 and count_hydrogens == 1:
            secondary_candidate += 1

    # If no candidate –OH group (non–acidic) was found, we cannot classify it.
    if total_candidate == 0:
        return False, "No identifiable non–acidic –OH groups found"
    
    if total_candidate == secondary_candidate:
        return True, f"All {total_candidate} candidate –OH group(s) are secondary (each on an sp³ carbon with 2 carbons and 1 hydrogen)"
    else:
        return False, (f"Not all identified –OH groups are secondary "
                       f"(found {secondary_candidate} secondary out of {total_candidate} candidate –OH groups)")

# Example usage:
# Uncomment these lines to run a couple of tests:
# test_list = [
#     "CCCCCCCCCCCCCCC(O)CCCC",  # nonadecan-5-ol: expected secondary alcohol.
#     "CC(O)CC",                # butan-2-ol: a textbook secondary alcohol.
#     "OC(=O)C(C/C=C\\C[C@H](\\C=C\\C=C/C=C/C=C/[C@H](C/C=C\\CC)O)O)O"  # aspirin-triggered resolvin D2 (acid OH should be ignored)
# ]
# for smi in test_list:
#     result, reason = is_secondary_alcohol(smi)
#     print(smi, "->", result, ":", reason)