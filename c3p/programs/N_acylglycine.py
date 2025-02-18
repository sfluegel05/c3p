"""
Classifies: CHEBI:16180 N-acylglycine
"""
"""
Classifies: N‐acylglycine 
Definition: An N‐acyl‐amino acid in which the amino acid is glycine.
This function attempts to verify that the molecule contains (at least one)
the substructure:
    R‑C(=O)‑N‑CH₂‑C(=O)O 
where the “CH₂” is enforced (by H2) to represent glycine.
In order to reduce false positives from peptides (which have multiple amide bonds),
we also count the total number of amide bonds ([NX3][CX3](=O)) in the molecule.
If more than one such amide bond is found, we assume that the match is part of
a peptide chain rather than an isolated N‐acylglycine.
"""

from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.

    An N-acylglycine is an N-acyl amino acid in which the amino acid is glycine.
    The canonical substructure is: 
          R-C(=O)-N-CH2-C(=O)O
    In this implementation we:
      1. Require that the CH2 group is exactly a methylene (i.e. has two hydrogens, H2).
      2. Use a SMARTS pattern that enforces the above.
      3. Count the number of amide bonds ([NX3][CX3](=O)); if more than one is found,
         we assume the molecule is a peptide (or has extra amide bonds) and thus not
         a simple N-acylglycine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an N-acylglycine, False otherwise.
        str: Explanation of the decision.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a refined SMARTS pattern for the N-acylglycine fragment.
    # [CX3](=O)   --> a carbonyl carbon (acyl carbon)
    # [NX3]       --> the amide nitrogen (trivalent)
    # [CH2;H2]    --> a methylene group (CH2 with exactly 2 hydrogens) representing glycine
    # [C](=O)[O]  --> a carboxyl group (COOH or COO-) attached to the glycine alpha carbon.
    #
    # This pattern does not enforce attachment to the remainder R (the acyl chain) or 
    # the acid proton; it just catches the overall connectivity.
    n_acylglycine_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3][CH2;H2][C](=O)[O]")
    if n_acylglycine_pattern is None:
        return False, "Failed to create SMARTS pattern for N-acylglycine"
    
    # Look for substructure matches of the N-acylglycine fragment.
    matches = mol.GetSubstructMatches(n_acylglycine_pattern)
    if not matches:
        return False, "N-acylglycine substructure not found"
    
    # To reduce false positives from peptides we count the number of amide bonds.
    # We define an amide bond SMARTS as any [NX3][CX3](=O) fragment.
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    n_amide_bonds = len(amide_matches)
    
    # In a simple N-acylglycine there should be only one amide bond.
    # If the molecule has additional amide bonds, it is likely a peptide or more complex molecule.
    if n_amide_bonds > 1:
        return False, f"Multiple amide bonds found ({n_amide_bonds}); likely part of a peptide chain"
    
    # If we have at least one match for the refined fragment and only one amide, we classify as N-acylglycine.
    return True, "Molecule contains the N-acylglycine substructure (R-C(=O)-N-CH2-C(=O)O)"

# Example usage (for testing independently):
# test_smiles = "CC(=O)NCC(O)=O"  # N-acetylglycine; expected True
# result, reason = is_N_acylglycine(test_smiles)
# print(result, reason)