"""
Classifies: CHEBI:33704 alpha-amino acid
"""
#!/usr/bin/env python
"""
Classifies: Alpha-amino acid
Definition: An amino acid in which the amino group is located on the carbon atom at the position alpha to the carboxy group.
"""

from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has an amino group attached to the alpha carbon, which also bears a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the alpha-amino acid motif.
    # Pattern 1: neutral carboxylic acid: looks for an amino group (with 1 or 2 hydrogens) attached to a carbon that in turn is attached to a -C(=O)O group.
    alpha_pattern_neutral = Chem.MolFromSmarts("[NX3;H2,H1][C;H](C(=O)O)")
    # Pattern 2: ionized carboxylate: same as above but carboxyl group depicted as C(=O)[O-].
    alpha_pattern_ionized = Chem.MolFromSmarts("[NX3;H2,H1][C;H](C(=O)[O-])")
    
    # Check if the molecule contains either pattern.
    if mol.HasSubstructMatch(alpha_pattern_neutral) or mol.HasSubstructMatch(alpha_pattern_ionized):
        return True, "Matches alpha-amino acid motif: amino group attached to the alpha carbon bearing a carboxyl (acid or carboxylate) group"
    else:
        return False, "Does not match the required pattern for an alpha-amino acid"

# Example usage (uncomment to run tests):
# test_smiles_list = [
#     "[H][C@@](\\N=C1CC(O)(CO)CC(NCC(O)=O)=C\\1OC)([C@H](C)O)C(O)=O",  # porphyra-334 (alpha-amino acid motif exists)
#     "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(O)=O",            # L-thyroxine (contains alpha-amino acid motif)
#     "CSCCCCCCC(N)C(O)=O",                                               # tetrahomomethionine (simple alpha-amino acid)
# ]
# for smi in test_smiles_list:
#     result, reason = is_alpha_amino_acid(smi)
#     print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")