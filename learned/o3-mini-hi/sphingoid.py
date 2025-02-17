"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: Sphingoid compounds
Definition: Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.
Improved strategy:
  - Split the molecule into connected fragments and require one fragment to satisfy:
    * A sphingoid core (one of several SMARTS variants) found entirely on acyclic atoms.
    * A long aliphatic chain (>=8 contiguous carbon atoms, with bonds that allow both single and double bonds) found by a relaxed SMARTS.
    * At least one nitrogen atom.
    * A fragment mass (roughly 200–800 Da) and carbon count (roughly 16–50).
  
False positives were arising because the global molecule could match the two criteria on two separate parts. Now we only classify a molecule if one fragment contains both features.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid compound based on its SMILES string.
    Sphingoid compounds are defined as sphinganine (and its homologs/stereoisomers) and their 
    hydroxy/unsaturated derivatives. For classification, we require that at least one connected 
    fragment (a candidate sphingoid fragment) satisfies:
      - It contains a sufficiently long aliphatic chain (allowing both saturated and unsaturated bonds)
        defined as at least 8 contiguous non‐ring carbons.
      - It contains at least one nitrogen (for the amino group).
      - It has at least one sphingoid core match (one of a few SMARTS patterns for a primary amino-diol 
        or its deoxy/carbonyl variants) whose atoms are all acyclic.
      - Its molecular weight is roughly between 200 and 800 Da and it has a reasonable total carbon count.
      
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        (bool, str): A tuple containing the classification (True if sphingoid, False otherwise) and an explanation.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Split molecule into disconnected fragments.
    frags = Chem.GetMolFrags(mol, asMols=True)
    if not frags:
        return False, "No fragments found in molecule."

    # Define a relaxed SMARTS pattern for a long aliphatic chain.
    # Matches a path of 8 contiguous carbons (allowing either sp3 or sp2) that are not in rings.
    chain_smarts = "[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]~[#6;!R]"
    chain_pattern = Chem.MolFromSmarts(chain_smarts)
    
    # Define a list of sphingoid core SMARTS patterns.
    core_smarts_list = [
        "C(O)C(N)CO",    # fully hydroxylated amino-diol core
        "C(=O)C(N)CO",   # carbonyl variant (dehydro form)
        "C(=O)CN",       # deoxy-carbonyl variant (shorter motif)
        "C(O)CN"         # deoxy-hydroxy variant
    ]
    core_patterns = [Chem.MolFromSmarts(s) for s in core_smarts_list if Chem.MolFromSmarts(s) is not None]

    # Helper: Check if the substructure match is entirely acyclic.
    def acyclic_match(mol_obj, pattern):
        matches = mol_obj.GetSubstructMatches(pattern)
        for match in matches:
            if all(not mol_obj.GetAtomWithIdx(idx).IsInRing() for idx in match):
                return True
        return False

    # Evaluate each fragment candidate.
    for frag in frags:
        # Compute fragment properties.
        frag_wt = rdMolDescriptors.CalcExactMolWt(frag)
        # Typical sphingoid compounds lie roughly between 200 and 800 Da.
        if frag_wt < 200 or frag_wt > 800:
            continue

        # Count carbons in fragment.
        c_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
        # Typical sphingoid fragments have a moderate-to-high carbon count.
        if c_count < 16 or c_count > 50:
            continue

        # Require at least one nitrogen in this fragment.
        n_count = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 7)
        if n_count < 1:
            continue

        # Check if the fragment has a long aliphatic chain.
        if not frag.HasSubstructMatch(chain_pattern):
            continue

        # Check if any of the sphingoid core patterns are found in an acyclic manner.
        core_found = any(acyclic_match(frag, core) for core in core_patterns)
        if not core_found:
            continue

        # If all tests passed for this fragment, classify as a sphingoid compound.
        return True, ("Fragment with weight {:.1f} Da, {} carbons, and {} nitrogen(s) "
                      "contains an acyclic sphingoid core and a long aliphatic chain."
                      .format(frag_wt, c_count, n_count))
        
    # If none of the fragments satisfied the conditions, return False.
    return False, "No fragment met the criteria for a sphingoid compound."

# For testing purposes, a few examples can be run:
if __name__ == "__main__":
    test_smiles = [
        # True positives (should be classified as sphingoid)
        "CCCCCCCCCCCC\\C=C\\C(=O)[C@@H](N)CO",   # 3-dehydrosphingosine
        "OC[C@@]([C@@](CCCCCCCCCCCCCC)(O)H)(N)H",  # Heptadecasphinganine
        # False negative previously (should now be accepted)
        "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CC[C@@H](O)[C@H](C)N",  # 1-deoxysphinga-6Z,9Z,12Z,15Z-tetraenine
        # A molecule that should not qualify.
        "CCO"
    ]
    for smi in test_smiles:
        result, reason = is_sphingoid(smi)
        print("SMILES:", smi)
        print("Is sphingoid?:", result)
        print("Reason:", reason)
        print()