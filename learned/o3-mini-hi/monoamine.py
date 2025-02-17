"""
Classifies: CHEBI:63534 monoamine
"""
"""
Classifies: Monoamine
Definition: A monoamine is an aralkylamino compound having exactly one amino group
connected to an aromatic ring by a two‐carbon chain. In this implementation we search 
for a specific SMARTS fragment that corresponds to a nitrogen (non‐aromatic, acyclic, and 
not adjacent to a carbonyl), bonded to two consecutive aliphatic carbons (non‐aromatic, non‐ring), 
which in turn is bonded to an aromatic carbon.
"""

from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.

    A monoamine is defined as an aralkylamino compound having exactly one amino group connected 
    to an aromatic ring by a two‐carbon chain. In our implementation we search for a substructure
    matching the SMARTS pattern:
      [#7;!a;!R;!$([#7]-C(=O))]-[#6;!a;!R]-[#6;!a;!R]-[c]
    which represents a nitrogen (that is non-aromatic, not in a ring, and not directly bonded to a carbonyl)
    bonded to an aliphatic carbon (first methylene), then to another aliphatic carbon (second methylene),
    and finally to an aromatic carbon.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule contains exactly one aralkylamino moiety as defined, False otherwise.
        str: A message explaining the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the aralkylamino moiety:
    #   [#7;!a;!R;!$([#7]-C(=O))] : a nitrogen atom that is not aromatic, not in a ring, and not adjacent to a carbonyl
    #   -[#6;!a;!R]           : an aliphatic (sp3) carbon, not aromatic and not in a ring
    #   -[#6;!a;!R]           : a second aliphatic carbon
    #   -[c]                  : an aromatic carbon
    smarts = "[#7;!a;!R;!$([#7]-C(=O))]-[#6;!a;!R]-[#6;!a;!R]-[c]"
    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        return False, "Error creating SMARTS pattern"

    # Find all matches of the pattern in the molecule.
    # Each match is a tuple of atom indices corresponding to the ordered atoms in the SMARTS.
    matches = mol.GetSubstructMatches(pattern)
    
    # To avoid counting duplicate discoveries of the same moiety (if more than one matching path walks over the same N)
    # we collect the unique nitrogen atom indices from the matches.
    unique_nitrogens = set()
    for match in matches:
        # match[0] is the index of the nitrogen in the matched fragment.
        unique_nitrogens.add(match[0])
        
    count = len(unique_nitrogens)
    
    # Verify exactly one aralkylamino moiety was found.
    if count == 0:
        return False, "No aralkylamino moiety (amine connected by a two‐carbon chain to an aromatic ring) found"
    elif count > 1:
        return False, f"Found {count} aralkylamino moieties; expected exactly one monoamine group"
    else:
        return True, "Contains exactly one amino group connected via a two‐carbon chain to an aromatic ring"

# Example usage (for manual testing):
if __name__ == "__main__":
    test_smiles = [
        "CNC[C@@H](O)c1ccc(O)c(O)c1",  # (S)-adrenaline: valid
        "NCCc1ccc(O)cc1",              # tyramine: valid
        "COc1cc(ccc1OCCCN1CCC(CC1)c1noc2cc(F)ccc12)C(C)=O",  # iloperidone: should be rejected
        "C(CNCCCCCCNCCC1=CC=CC=C1)C2=CC(O)=C(C=C2)O"  # dopexamine: has two aralkylamino moieties; should be rejected
    ]
    for s in test_smiles:
        result, reason = is_monoamine(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")