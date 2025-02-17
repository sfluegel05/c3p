"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: prenylquinone

A prenylquinone is defined as a quinone substituted by a polyprenyl‐derived side chain.
That is, the molecule must contain a quinone core (for example, a 1,4‐benzoquinone or a naphthoquinone motif)
and a prenyl-derived side chain that contains at least two isoprene units (e.g. two or more “[CH2]C=C([CH3])” fragments)
directly attached to the quinone core.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    
    The molecule must contain:
      1. A quinone core represented by one of several common SMARTS patterns.
      2. At least one substituent branch on that core that is prenyl-derived,
         and that branch must have at least two isoprene units (heuristic: two or more matches
         of the motif [CH2]C=C([CH3]) that are directly attached to the core).
      
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule is classified as a prenylquinone, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1: Find the quinone core --- 
    # We use a list of SMARTS patterns that cover several quinone types (benzoquinone, naphthoquinone, etc.)
    quinone_patterns = [
        "c1cc(=O)cc(=O)c1",                # 1,4-benzoquinone
        "c1ccc2C(=O)c(c1)C(=O)cc2",         # naphthoquinone variant 1
        "c1cc2c(c(c1)C(=O))C(=O)cc2"         # naphthoquinone variant 2
    ]
    quinone_core_indices = set()
    for qs in quinone_patterns:
        pattern = Chem.MolFromSmarts(qs)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            # For our purposes, we take the atoms in the first matching pattern as the quinone core.
            quinone_core_indices = set(matches[0])
            break
            
    if not quinone_core_indices:
        return False, "No quinone core detected; none of the expected quinone SMARTS matched"
    
    # --- Step 2: Identify prenyl fragment(s) ---
    # Use a SMARTS pattern to capture an isoprene unit.
    prenyl_pattern = Chem.MolFromSmarts("[CH2]C=C([CH3])")
    if prenyl_pattern is None:
        return False, "Error in prenyl SMARTS pattern"
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)
    if not prenyl_matches:
        return False, "No prenyl fragment detected"
    
    # --- Step 3: Check if at least one prenyl-derived substituent (with at least 2 isoprene units)
    # is directly attached to the quinone core. We do this by:
    #   (a) Identifying prenyl matches that are directly bonded (at least one atom) to an atom in the quinone core.
    #   (b) Counting, per substituent branch, how many distinct prenyl matches occur.
    attached_prenyl_counts = 0
    for match in prenyl_matches:
        # For each prenyl match, check if any atom is adjacent to a quinone core atom.
        is_attached = False
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in quinone_core_indices:
                    is_attached = True
                    break
            if is_attached:
                break
        if is_attached:
            attached_prenyl_counts += 1
    
    # For a true prenylquinone the polyprenyl side chain should have at least two i.e. 2 or more isoprene units.
    if attached_prenyl_counts < 2:
        return False, f"Prenyl fragment detected but only {attached_prenyl_counts} isoprene unit(s) are attached to the quinone core (need at least 2 for a polyprenyl side chain)"
    
    return True, f"Molecule contains a quinone core with {attached_prenyl_counts} prenyl (isoprene) units attached"

# Example test run
if __name__ == "__main__":
    # Use ubiquinone-2 as one of the examples.
    test_smiles = "COC1=C(OC)C(=O)C(C\\C=C(/C)CCC=C(C)C)=C(C)C1=O"
    result, reason = is_prenylquinone(test_smiles)
    print(f"Result: {result}, Reason: {reason}")