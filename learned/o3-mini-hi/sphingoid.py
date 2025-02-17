"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: Sphingoid compounds
Definition: Sphinganine, its homologs and stereoisomers, and the hydroxy and unsaturated derivatives of these compounds.
This version improves core detection by considering several SMARTS patterns for an acyclic sphingoid core.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid compound based on its SMILES string.
    Sphingoid compounds are defined as sphinganine (and its homologs/stereoisomers) and their
    hydroxy/unsaturated derivatives. In our classifier, the molecule is expected to have:
      - A long aliphatic (acyclic) chain of at least 8 contiguous sp3 carbons.
      - At least one nitrogen atom (from the amino group).
      - A characteristic sphingoid core featuring an amino-diol (or deoxy/carbonyl variant) motif.
        Instead of a single motif we now check a set of similar patterns:
          Pattern 1: C(O)C(N)CO        (fully hydroxylated core)
          Pattern 2: C(=O)C(N)CO       (carbonyl variant)
          Pattern 3: C(=O)CN           (deoxy-carbonyl variant)
          Pattern 4: C(O)CN            (deoxy-hydroxy variant)
        At least one match must be found completely on an acyclic fragment.
      - A molecular weight above ~200 Da.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as sphingoid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a sufficient long alkyl chain.
    # We require a contiguous chain of at least 8 carbon atoms (SMARTS: "CCCCCCCC").
    chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No sufficiently long aliphatic chain (>=8 contiguous carbons) found."
    
    # Check total carbon count (to avoid tiny molecules).
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, f"Too few carbons ({c_count}); sphingoid molecules usually have long aliphatic chains."
    
    # Require at least one nitrogen (for the amino group).
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "No nitrogen atom found; sphingoid compounds require an amino group."
    
    # List a set of SMARTS patterns to capture variations of the sphingoid core.
    core_smarts_list = [
        "C(O)C(N)CO",   # fully hydroxylated amino-diol core
        "C(=O)C(N)CO",  # carbonyl variant (dehydro form)
        "C(=O)CN",      # deoxy-carbonyl variant (shorter motif)
        "C(O)CN"        # deoxy-hydroxy variant
    ]
    core_patterns = [Chem.MolFromSmarts(smarts) for smarts in core_smarts_list]
    
    # Function to check if any match of a given pattern is found on an entirely acyclic fragment.
    def acyclic_match(pattern):
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            if all(not mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
                return True
        return False

    core_found = any(acyclic_match(core) for core in core_patterns if core is not None)
    if not core_found:
        return False, "No acyclic sphingoid core (amino-diol or its deoxy/carbonyl variant) found."
    
    # Check molecular weight: Most sphingoid molecules lie roughly between 200 and 800 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a typical sphingoid compound."
    
    # All tests passed: classify as sphingoid.
    return True, "Molecule exhibits an acyclic sphingoid core with a long aliphatic chain and appropriate functionality."

# Example test cases (for manual testing)
if __name__ == "__main__":
    test_smiles = [
        "CCCCCCCCCCCC\\C=C\\C(=O)[C@@H](N)CO",  # 3-dehydrosphingosine (should return True)
        "OC[C@@]([C@@](CCCCCCCCCCCCCC)(O)H)(N)H",  # Heptadecasphinganine (True)
        "CCO",  # Too short (False)
        "CCCCCCCCCCCCCCC(=O)CN",  # 1-deoxymethyl-3-dehydrosphinganine (should now match Pattern 3)
    ]
    for smi in test_smiles:
        classified, explanation = is_sphingoid(smi)
        print("SMILES:", smi)
        print("Classified as sphingoid:", classified)
        print("Reason:", explanation)
        print()