"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N‐acylglycine
Definition: An N‐acyl‐amino acid in which the amino acid specified is glycine.
Canonical fragment: R-C(=O)-N-CH2-C(=O)[O]
Requirements:
  • There is an acyl carbonyl group attached to an amide nitrogen.
  • The amide nitrogen is attached to a glycine –CH2– group,
    which in turn is attached to a carboxyl carbon.
  • The amide nitrogen must only be attached to the acyl carbon and the glycine carbon.
  • The glycine carbon is a CH2 (exactly 2 hydrogens) and has no extra heavy substitutes.
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.

    An N-acylglycine is an N-acyl amino acid where the amino acid is glycine.
    Its canonical fragment is: R-C(=O)-N-CH2-C(=O)[O]
    (the terminal carboxyl group may appear as -COO⁻ or -COOH).

    Process:
      1. Parse the SMILES and add explicit hydrogens.
      2. Try two SMARTS patterns to capture both ionized and unionized forms:
             Pattern A: [CX3](=O)[NX3][CH2;H2][CX3](=O)[O-]
             Pattern B: [CX3](=O)[NX3][CH2;H2][CX3](=O)[OH]
      3. For each match, verify:
             a. The amide nitrogen (match atom 1) is bonded only to the acyl carbon and the glycine carbon.
             b. The glycine carbon (match atom 2) is a CH2 (exactly two heavy neighbors: one is the amide nitrogen, one is the carboxyl carbon)
                and has exactly 2 hydrogens (using GetTotalNumHs).
      4. If a valid match is found, return True with explanation; otherwise return False.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an N-acylglycine, False otherwise.
        str: Explanation of the decision.
    """
    # Parse the molecule from SMILES and add explicit hydrogens for accurate hydrogen count.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)

    # Define two SMARTS patterns to cover both deprotonated (carboxylate) and protonated carboxyl groups.
    smarts_list = [
        "[CX3](=O)[NX3][CH2;H2][CX3](=O)[O-]",  # deprotonated version
        "[CX3](=O)[NX3][CH2;H2][CX3](=O)[OH]"     # protonated version
    ]

    for smarts in smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue  # skip if pattern creation fails
        
        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            continue  # try next SMARTS
        # Each match is a tuple of indices corresponding to:
        # 0: Acyl carbon (carbonyl of the acyl group)
        # 1: Amide nitrogen
        # 2: Glycine methylene carbon (should be CH2)
        # 3: Carboxyl carbon (of the glycine acid group)
        for match in matches:
            acyl_c = mol.GetAtomWithIdx(match[0])
            amide_n = mol.GetAtomWithIdx(match[1])
            glycine_c = mol.GetAtomWithIdx(match[2])
            carboxyl_c = mol.GetAtomWithIdx(match[3])
            
            # Check connectivity for the amide nitrogen:
            # It should have exactly two heavy (non-H) neighbors: one acyl carbon and one glycine carbon.
            n_heavy = [nbr for nbr in amide_n.GetNeighbors() if nbr.GetAtomicNum() > 1]
            if set(n_heavy) != {acyl_c, glycine_c}:
                continue   # Not matching required connectivity
            
            # Check glycine carbon connectivity:
            # It should have exactly two heavy neighbors: the amide nitrogen and the carboxyl carbon.
            glycine_heavy = [nbr for nbr in glycine_c.GetNeighbors() if nbr.GetAtomicNum() > 1]
            if set(glycine_heavy) != {amide_n, carboxyl_c}:
                continue   # Extra connections on glycine carbon
            
            # Also ensure that the glycine carbon really is a CH2 (with exactly 2 hydrogens)
            if glycine_c.GetTotalNumHs() != 2:
                continue
            
            # If all checks pass, we consider the substructure as valid.
            return True, "Molecule contains a valid N-acylglycine fragment: R-C(=O)-N-CH2-C(=O)[O]"

    # If no match satisfies all connectivity conditions, then return failure.
    return False, "N-acylglycine substructure not found or key atom connectivity is ambiguous."

# Example usage (uncomment to test):
# test_smiles = "CC(=O)NCC(O)=O"  # Should be True for N-acetylglycine.
# result, reason = is_N_acylglycine(test_smiles)
# print(result, reason)