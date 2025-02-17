"""
Classifies: CHEBI:33704 alpha-amino acid
"""
#!/usr/bin/env python
"""
Classifies: Alpha–amino acid
Definition: An amino acid in which the amino group is located on the carbon atom 
at the position alpha to the carboxylic acid.
    
Our improved approach:
  • First, we parse the molecule (adding explicit H’s so that the N’s hydrogen count is reliable).
  • We reject molecules that contain any peptide (amide) bond (which would signal a di-/oligo–peptide).
  • Then we use a SMARTS pattern that attempts to capture a “naked” alpha–amino acid motif. 
    In a canonical amino acid the alpha carbon (Cα) is bound to:
      – exactly one hydrogen (so we require [C;H1]),
      – an amino group (we use [NX3] and then later check that at least one hydrogen is present on that N),
      – and a carboxyl group. For the carboxyl group we search for a carbon (C) doubly bonded to 
        one oxygen and singly bonded to an oxygen that carries a hydrogen (or is deprotonated).
    The SMARTS "[C;H1]([NX3])[C](=O)[O;H1,O-]" encodes that notion.
  • Finally, we inspect the match to ensure that (a) the carboxyl carbon is truly terminal – that is, 
    it must have only two oxygen neighbors (besides the connection to Cα) – and (b) the amino nitrogen 
    carries at least one hydrogen (to avoid quaternary N’s).
  • We also require that exactly one such motif is found in the molecule.
    
This approach is not bulletproof but should better separate the true alpha–amino acids from the false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid(smiles: str):
    """
    Determines whether a molecule (given by its SMILES) is a stand-alone alpha–amino acid.
    
    The algorithm:
      1. Parses and adds explicit hydrogens.
      2. Rejects molecules that contain peptide (amide) bonds.
      3. Uses a SMARTS query to search for the core motif:
             [C;H1]([NX3])[C](=O)[O;H1,O-]
         This means: an sp3 (or nearly tetrahedral) alpha-carbon carrying exactly one non–H, that is bound 
         to a nitrogen (whose hydrogen count will be verified) and to a carboxyl carbon that is bound to two oxygens.
      4. For the matched motif, do additional local checks:
            - The carboxyl carbon should have exactly two oxygen neighbors (the carbonyl and the hydroxy oxygen)
              (i.e. be “terminal”).
            - The amino nitrogen should have at least one implicit/explicit hydrogen.
      5. Accept only if exactly one such valid motif is found.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): Tuple of (classification, explanation).
    """
    # Parse the molecule and add explicit hydrogens for reliable hydrogen counting.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    molH = Chem.AddHs(mol)
    
    # Reject molecules that have peptide (amide) bonds.
    # The SMARTS "C(=O)N[C]" is a simple probe for a peptide bond.
    peptide_smarts = Chem.MolFromSmarts("C(=O)N[C]")
    if molH.HasSubstructMatch(peptide_smarts):
        return False, "Contains peptide bond(s) indicating a peptide rather than a single amino acid"
    
    # Define the SMARTS for the alpha–amino acid core motif.
    # It requires an alpha carbon with one hydrogen ([C;H1]),
    # which is bound to a nitrogen ([NX3]) and to a carboxyl carbon [C](=O)[O;H1,O-].
    aa_smarts = Chem.MolFromSmarts("[C;H1]([NX3])[C](=O)[O;H1,O-]")
    if aa_smarts is None:
        return None, None  # if SMARTS fails to compile

    matches = molH.GetSubstructMatches(aa_smarts)
    if len(matches) == 0:
        return False, "No alpha–amino acid motif found"
    if len(matches) > 1:
        return False, f"Multiple potential alpha–amino acid motifs found ({len(matches)}); not a simple amino acid"
    
    # Exactly one match was found.
    # The match is a tuple of atom indices. We assume:
    #  index 0: the alpha carbon,
    #  index 1: the attached nitrogen,
    #  index 2: the carboxyl carbon.
    match = matches[0]
    alpha_idx, n_idx, carboxyl_idx, *rest = match  # note: the SMARTS has 3 explicit atoms but RDKit may return more indices if the pattern is extended
    # Retrieve the atoms.
    alpha_atom = molH.GetAtomWithIdx(alpha_idx)
    n_atom = molH.GetAtomWithIdx(n_idx)
    carboxyl_atom = molH.GetAtomWithIdx(carboxyl_idx)
    
    # Verify that the carboxyl carbon is a terminal carboxyl.
    # Count oxygen neighbors (ignoring hydrogens).
    oxygens = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(oxygens) != 2:
        return False, "Carboxyl carbon does not appear to be terminal (should have exactly 2 oxygen neighbors)"
    
    # Verify that the amino nitrogen has at least one hydrogen (if not in a ring, we require at least one H)
    if n_atom.GetTotalNumHs() < 1:
        return False, "The amino nitrogen does not carry any hydrogen (not a free amino group)"
    
    # (Optional) We could require that the alpha carbon (ignoring its hydrogens) has exactly 3 non–H heavy neighbors:
    heavy_neighbors = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    # In a canonical amino acid the alpha carbon will be bound to: the amino nitrogen, the carboxyl carbon, and one side–chain atom.
    if len(heavy_neighbors) != 3:
        return False, f"Alpha carbon has {len(heavy_neighbors)} heavy atom neighbors rather than the expected 3"
    
    # Our “core” motif has been located and validated.
    return True, ("Matches alpha–amino acid motif: a free carboxyl group attached to an alpha carbon that "
                   "bears an amino group (with at least one hydrogen) and no peptide bonds are present")

    
# Example usage:
if __name__ == "__main__":
    # List a few examples from the provided list.
    test_smiles = [
        "[H][C@@](\\N=C1CC(O)(CO)CC(NCC(O)=O)=C\\1OC)([C@H](C)O)C(O)=O",  # porphyra-334 (should be true)
        "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(O)=O",            # L-thyroxine (should be true)
        "NC1(CC1)C(O)=O",                                                # 1-aminocyclopropanecarboxylic acid (true)
        "[O-]C(=O)[C@H]1CNCC[NH2+]1",                                     # (R)-piperazine-2-carboxylic acid zwitterion (false positive in previous attempt)
    ]
    for smi in test_smiles:
        result, reason = is_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")