"""
Classifies: CHEBI:26607 saturated fatty acid
"""
"""
Classifies: Saturated Fatty Acid
Definition: Any fatty acid containing no carbon–carbon multiple bonds.
They are characterized by a (terminal) carboxyl group and a fully saturated (i.e. no C=C or C≡C outside of the carbonyl)
hydrocarbon chain.
Additional guard conditions (allowed elements, no aromatic systems) are used to reduce false positives.
"""

from rdkit import Chem

def is_saturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a saturated fatty acid based on its SMILES string.
    
    In this approach we require:
      1. The molecule can be parsed.
      2. All heavy atoms are among allowed elements (here: C, H, O, and D for deuterium) – this helps weed out peptides etc.
      3. The molecule has no aromatic atoms.
      4. There is at least one carboxyl (or acyl) group. For this to be “fatty‐acid–like”
         the carboxyl carbon must be terminal (attached to exactly one carbon, its α–carbon).
      5. Outside of the C=O of the acid group, no other C=C or C≡C bonds exist.
         
    Note: Some lipids (e.g. fatty acyl chains in glycerophospholipids) may be missed because the acid group is masked 
    (or the molecule contains extra elements).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a saturated fatty acid, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Check that only allowed elements are present.
    # Allowed: Carbon, Hydrogen, Oxygen, and deuterium (“D”) (as a label on H).
    allowed = {"C", "H", "O", "D"}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in allowed:
            return False, f"Disallowed element found ({symbol}); not a simple fatty acid"
    
    # Step 2. If any atoms are aromatic, reject.
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic():
            return False, "Contains aromatic system; not a typical saturated fatty acid"
    
    # Step 3. Look for a carboxylic acid substructure.
    # SMARTS: C(=O)[O;H,-] matches a carbon with a double-bonded oxygen and an oxygen that is either neutral (with H) or anionic.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H,-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # For a fatty acid the acid carbon (first atom in pattern) should be terminal;
    # i.e. it should have exactly one carbon neighbor (its α–carbon).
    terminal_acid_found = False
    for match in acid_matches:
        acid_carbon_idx = match[0]
        acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
        # Count number of carbon neighbors attached (ignore the double-bonded O or the -OH).
        carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_acid_found = True
            break
    if not terminal_acid_found:
        return False, "Carboxylic acid group is not terminal (acid carbon has >1 C neighbor)"
    
    # Step 4. Check for other carbon–carbon multiple bonds.
    # We do not want any C=C or C≡C bonds outside of the allowed acid carbonyl.
    # Use a SMARTS that matches a C=C bond when neither carbon is attached to an oxygen by a double bond.
    unsat_double = Chem.MolFromSmarts("[$([C;!$(*=O)]):1]=[$([C;!$(*=O)]):2]")
    if mol.HasSubstructMatch(unsat_double):
        return False, "Contains carbon–carbon double bond(s) outside of carboxyl group"
    unsat_triple = Chem.MolFromSmarts("[$([C;!$(*=O)]):1]#[C]")
    if mol.HasSubstructMatch(unsat_triple):
        return False, "Contains carbon–carbon triple bond(s)"
    
    # (Optionally, one could add a check that the molecule is not extremely “small” or that the ratio of C:O is as expected 
    # for a fatty acid; however, for example propionic acid (CCC(O)=O) is a valid saturated fatty acid.)
    
    return True, "Saturated fatty acid: contains terminal carboxylic acid group and no carbon–carbon multiple bonds"

# (Optional main/test block could be added here to try some SMILES.)
if __name__ == "__main__":
    test_smiles = [
        "CC(=O)CCC(O)=O",  # 4-oxopentanoic acid, true positive
        "CCCC(C)C(O)=O",   # 2-methylvaleric acid, true positive
        "OC(=O)CCC(CCCC)CC",  # 4-Ethyloctanoic acid (false positive in previous code)
        "C(C[N+](C)(C)C)OP(=O)([O-])OC[C@H](OC)COC(CCCCCCCCCCCCCCC)=O"  # phosphocholine type (false negative previously)
    ]
    for s in test_smiles:
        result, reason = is_saturated_fatty_acid(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n{'-'*40}")