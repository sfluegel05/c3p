"""
Classifies: CHEBI:33704 alpha-amino acid
"""
#!/usr/bin/env python
"""
Classifies: Alpha-amino acid
Definition: An amino acid in which the amino group is located on the carbon atom at the position alpha to the carboxy group.
We improve on the previous approach by:
  • Rejecting molecules that contain peptide (amide) bonds.
  • Searching for an “alpha carbon” that is directly linked to any nitrogen atom and a free acid or carboxylate group.
This heuristic is not perfect but helps reduce false positives (e.g. peptides) and also catch valid amino acids that lack explicit chirality.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    In this revised approach the molecule is accepted as a stand-alone alpha-amino acid if:
      1. It does not contain any peptide bonds (i.e. no C(=O)N substructure),
         which helps weed out peptides with free amino termini.
      2. It contains at least one occurrence of an “alpha amino acid” motif defined as
         a tetrahedral carbon attached (directly) both to at least one nitrogen (any N)
         and to a carboxyl group (either shown as C(=O)O or C(=O)[O-]). 
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a (stand-alone) alpha-amino acid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules that contain peptide bonds.
    # (This catches dipeptides, tripeptides, etc. even if an N-terminus might match the motif.)
    amide_bond_smarts = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_bond_smarts):
        return False, "Contains peptide bond(s) indicating a peptide rather than a single amino acid"

    # Define SMARTS for free carboxyl groups.
    free_carboxyl_patterns = [
        Chem.MolFromSmarts("C(=O)O"),    # neutral acid form
        Chem.MolFromSmarts("C(=O)[O-]")   # deprotonated carboxylate form
    ]
    
    # Define SMARTS for the alpha-amino acid motif.
    # The central alpha carbon should be attached to a nitrogen (in any oxidation state) and a carboxyl group.
    # We do not specify chirality or the number of attached hydrogens so that variations (e.g. N-methyl, ring systems)
    # are accepted.
    alpha_patterns = []
    for fc in free_carboxyl_patterns:
        # Pattern: [C]([#7])(C(=O)O or C(=O)[O-])
        # Using the wildcard atomic number for nitrogen: [#7]
        # The [C] is unconstrained aside from being the attachment point.
        patt = Chem.MolFromSmarts(f"[C]([#7])({Chem.MolToSmarts(fc)})")
        if patt is not None:
            alpha_patterns.append(patt)
    
    # Search for at least one match to the alpha amino acid motif.
    for patt in alpha_patterns:
        if mol.HasSubstructMatch(patt):
            return True, "Matches alpha-amino acid motif: a carbon with a bonded nitrogen and a free carboxyl group, with no peptide bonds present"
    
    return False, "Does not match the required alpha-amino acid motif"

# Example test cases (uncomment to run some tests):
# test_smiles = [
#     "[H][C@@](\\N=C1CC(O)(CO)CC(NCC(O)=O)=C\\1OC)([C@H](C)O)C(O)=O",  # porphyra-334, should be True
#     "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(O)=O",            # L-thyroxine, should be True
#     "NC1(CC1)C(O)=O",                                               # 1-aminocyclopropanecarboxylic acid, should be True
#     "CNCC(O)=O",                                                   # sarcosine, should be True
#     "O=C(NCC(=O)N[C@@H](CO)C(O)=O)[C@@H](N)CC=1NC=NC1",             # His-Gly-Ser (dipeptide/peptide), should be False
#     "S(CC[C@@](N)(C)C(O)=O)C"                                        # alpha-methylmethionine, should be True
# ]
#
# for smi in test_smiles:
#     result, reason = is_alpha_amino_acid(smi)
#     print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")