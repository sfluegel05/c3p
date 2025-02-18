"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: glycerophosphocholine
Definition: The glycerol phosphate ester of a phosphocholine.
A true glycerophosphocholine should contain (i) a phosphocholine headgroup – i.e.
an oxygen-bridged phosphate with a two-carbon chain ending in a trimethylammonium group –
and (ii) that phosphate must be connected via an oxygen (typically on a primary CH2-bearing carbon)
to a short glycerol (or glycerol-like) fragment. In many phospholipids (e.g. diacyl phosphatidylcholines)
the glycerol bears two acyl chains, but those are considered false positives here.
This heuristic method uses SMARTS patterns that capture (a) the ideal non-acylated glycerophosphocholine
and (b) variants with exactly one acyl ester substituent on the glycerol.
Note: This method is heuristic and may not perfectly capture all valid (or invalid) cases.
"""

from rdkit import Chem

def is_glycerophosphocholine(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is a glycerophosphocholine.
    
    The approach is to look for a substructure that unites the required elements:
      (a) a phosphocholine headgroup: an oxygen-linked phosphate that is further connected to a 
          two-carbon chain terminating with a trimethylammonium,
      (b) connectivity to a glycerol or glycerol-like backbone—that is, a three-carbon fragment where
          the phosphate is esterified to one of the primary (CH2) groups, and at most one acyl (ester) group 
          is present on the glycerol.
      
    We define three SMARTS patterns:
      1. A fully “free” glycerol backbone (i.e. no acyl chain).
         The SMARTS below expresses the connectivity:
           HO-CH2–CHOH–CH2–O–P(=O)(O)–O–CH2CH2–N⁺(C)(C)C
         as: "OCC(O)CO[P](=O)(O)OCC[N+](C)(C)C"
      2. A lysophosphatidylcholine pattern with an acyl at the sn-1 position:
           Acyl at sn-1 means the leftmost hydroxyl (i.e. HO– becomes O–C(=O)...)
         Pattern: "OC(=O)OCC(O)CO[P](=O)(O)OCC[N+](C)(C)C"
      3. A lysophosphatidylcholine pattern with an acyl at the sn-2 position:
           Acyl at sn-2 means the central hydroxyl is acylated.
         Pattern: "OCC(OC(=O)*)CO[P](=O)(O)OCC[N+](C)(C)C"
    
    If any of these patterns are found, we accept the molecule as a glycerophosphocholine.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is accepted as glycerophosphocholine, otherwise False.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # List of SMARTS patterns to match the glycerophosphocholine connectivity.
    # Pattern 1: fully free glycerol backbone.
    pattern_free = "OCC(O)CO[P](=O)(O)OCC[N+](C)(C)C"
    # Pattern 2: lysophosphatidylcholine with one acyl chain at position sn-1.
    pattern_lyso_sn1 = "OC(=O)OCC(O)CO[P](=O)(O)OCC[N+](C)(C)C"
    # Pattern 3: lysophosphatidylcholine with one acyl chain at position sn-2.
    # In the SMARTS below, the '*' after OC(=O) allows for any carbon chain.
    pattern_lyso_sn2 = "OCC(OC(=O)*)CO[P](=O)(O)OCC[N+](C)(C)C"
    
    patterns = [
        (pattern_free, "Found free glycerol backbone connectivity"),
        (pattern_lyso_sn1, "Found lysophosphatidylcholine connectivity (acyl at sn-1)"),
        (pattern_lyso_sn2, "Found lysophosphatidylcholine connectivity (acyl at sn-2)")
    ]
    
    for patt_str, explanation in patterns:
        patt = Chem.MolFromSmarts(patt_str)
        if patt is None:
            continue
        if mol.HasSubstructMatch(patt):
            return True, explanation
    
    return False, "Glycerophosphocholine connectivity not found; either phosphocholine headgroup or proper glycerol linkage is missing or has too many acyl substitutions."

# Example usage:
if __name__ == "__main__":
    # Test with one example provided.
    test_smiles = "P(OC[C@@H](CO)OC(CCCCCCC/C=C\\CC(CCCCCC)O)=O)(=O)(OCC[N+](C)(C)C)[O-]"
    result, reason = is_glycerophosphocholine(test_smiles)
    print(result, reason)