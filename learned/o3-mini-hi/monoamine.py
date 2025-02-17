"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: Monoamine, defined as an aralkylamino compound which contains one amino group 
connected to an aromatic ring by a two‐carbon chain.
Examples include tyramine, (R)-noradrenaline, dopamine, etc.
"""

from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is defined as an aralkylamino compound having one amino group connected
    to an aromatic ring by a two-carbon chain. This means that there must be a substructure
    where a non-aromatic nitrogen (of any substitution) is bonded to two non-aromatic carbons,
    the second of which is bonded to an aromatic carbon.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for an aralkylamino group:
    # "[NX3;!a]" : a non-aromatic (aliphatic) trivalent nitrogen (could be primary, secondary or tertiary)
    # "[#6;!a]"  : a non-aromatic carbon atom (any carbon in a two-carbon chain)
    # "[c]"      : an aromatic carbon atom.
    # Thus the pattern "[NX3;!a][#6;!a][#6;!a][c]" finds a fragment N--C--C--(Ar).
    smarts = "[NX3;!a][#6;!a][#6;!a][c]"
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        return False, "Error in SMARTS pattern"
    
    # Find all substructure matches that satisfy the pattern.
    matches = mol.GetSubstructMatches(pattern)
    
    # Check if exactly one aralkylamino moiety was found.
    if len(matches) == 0:
        return False, "No aralkylamino moiety (amine connected by a two‐carbon chain to an aromatic ring) found"
    elif len(matches) > 1:
        return False, f"Found {len(matches)} aralkylamino moieties; expected exactly one monoamine group"
  
    # If exactly one match is found, we classify the molecule as a monoamine.
    return True, "Contains exactly one amino group connected via a two-carbon chain to an aromatic ring"

# Example test cases (uncomment to run)
# test_smiles = [
#     "NCCc1ccc(O)cc1",  # tyramine
#     "CNC[C@@H](O)c1ccc(O)c(O)c1",  # (R)-adrenaline
#     "C1=CC(O)=C(C=C1)CCN",  # simplified example that should match
#     "CC(C)NC[C@H](O)c1ccc(O)c(O)c1"  # L-isoprenaline
# ]
# for s in test_smiles:
#     result, reason = is_monoamine(s)
#     print(f"SMILES: {s} -> {result} ({reason})")