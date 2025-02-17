"""
Classifies: CHEBI:33704 alpha-amino acid
"""
#!/usr/bin/env python
"""
Classifies: Alpha–amino acid
Definition: An amino acid in which the amino group is located on the carbon atom at the position alpha to the carboxy group.
    
This implementation improves on the last try by:
  • Adding explicit hydrogens to count atoms more reliably.
  • Rejecting molecules that contain peptide (amide) bonds.
  • Using a relaxed SMARTS pattern that requires an amino group (any nitrogen) bound to a carbon (the α–carbon)
    which is in turn bound to a carboxyl group (C(=O)[O] where the oxygen can be in its protonated or deprotonated form).
  • Post–processing the motif to require that the carboxyl carbon is “terminal” (i.e. has exactly two oxygen heavy neighbors)
    and that the α–carbon has 2–3 heavy neighbors (the minimal connectivity expected for a stand–alone amino acid).
    
Note: This approach will not capture every possible edge case but should do a better job on the provided examples.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string) is a stand-alone alpha–amino acid.
    
    The algorithm:
      1. Parse the molecule and add explicit hydrogens.
      2. Reject the molecule if any peptide (amide) bonds are detected.
      3. Search for the core motif by using this SMARTS:
             [#7][C]([*])C(=O)[O;H1,O-]
         This pattern means: a nitrogen (any, with no explicit hydrogen constraint) attached to a carbon (our α–carbon)
         that in turn is bound to a carboxyl carbon (C(=O)[O] where the “hydroxy” oxygen can be either protonated or negative).
         We add a dummy neighbor ([*]) on the α–carbon to force it to be “substituted” by a side–chain (or hydrogen in glycine).
      4. For each match (if any), verify:
            - That exactly one motif is present (if more than one, then the molecule is likely not a single amino acid).
            - That the carboxyl carbon is terminal – i.e. it has exactly 2 oxygen neighbors.
            - That the α–carbon’s heavy–atom connectivity is acceptable: for glycine it is 2 (N and COO) and for others normally 3.
      5. If a valid motif is found, we return True with an explanation; otherwise, return False along with the reason.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): Tuple (True, explanation) if the molecule is a (stand–alone) alpha–amino acid,
                   otherwise (False, explanation). If the task cannot be performed, returns (None, None).
    """
    # Parse the molecule and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    molH = Chem.AddHs(mol)
    
    # Reject molecules that contain peptide (amide) bonds.
    # We look for the pattern "C(=O)N[C]" to flag amide bonds linking amino acids.
    peptide_smarts = Chem.MolFromSmarts("C(=O)N[C]")
    if molH.HasSubstructMatch(peptide_smarts):
        return False, "Contains peptide bond(s) indicating a peptide rather than a single amino acid"
    
    # Define the SMARTS for the alpha–amino acid core motif.
    # This pattern looks for:
    #   1) A nitrogen atom ([#7]) that can be bound as a free or substituted amine.
    #   2) An α–carbon ([C]) attached (via any neighbor, [*]) to enforce some substitution.
    #   3) Followed by a carboxyl group: a carbon with a double bond to oxygen and a single bond to an -OH (or O–).
    # Note that this pattern has four atoms: (N) - (α–C) - (COO–C) - (one oxygen of the carboxyl).
    aa_smarts = Chem.MolFromSmarts("[#7][C]([*])C(=O)[O;H1,O-]")
    if aa_smarts is None:
        return None, None  # SMARTS could not compile
    
    matches = molH.GetSubstructMatches(aa_smarts)
    if not matches:
        return False, "No alpha–amino acid motif found"
    if len(matches) > 1:
        return False, f"Multiple potential alpha–amino acid motifs found ({len(matches)}); not a simple amino acid"
    
    # Exactly one match was found.
    # The order of atom indices in our SMARTS is assumed to be:
    #   0: Amino nitrogen,
    #   1: α–carbon,
    #   2: Carboxyl carbon,
    #   3: One oxygen of the carboxyl group.
    match = matches[0]
    if len(match) < 3:
        return False, "Matched motif does not contain enough atoms"
    amino_idx = match[0]
    alpha_idx = match[1]
    carboxyl_idx = match[2]
    
    # Retrieve the atoms.
    amino_atom = molH.GetAtomWithIdx(amino_idx)
    alpha_atom = molH.GetAtomWithIdx(alpha_idx)
    carboxyl_atom = molH.GetAtomWithIdx(carboxyl_idx)
    
    # Check that the carboxyl carbon appears terminal. It should have exactly two oxygen neighbors
    # (one for the carbonyl and one for the hydroxy – which may be deprotonated).
    oxy_neighbors = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
    if len(oxy_neighbors) != 2:
        return False, "Carboxyl carbon does not appear to be terminal (should have exactly 2 oxygen neighbors)"
    
    # For the α–carbon, check its heavy atom (non–hydrogen) connectivity.
    # In a stand–alone amino acid, the α–carbon is attached to the amino group and the carboxyl group,
    # plus one side–chain (which in the case of glycine may be another hydrogen instead).
    heavy_neighbors = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    if len(heavy_neighbors) not in (2, 3):
        return False, f"Alpha carbon has {len(heavy_neighbors)} heavy atom neighbors instead of the expected 2 or 3"
    
    # If we reached this point, we believe the core alpha–amino acid motif has been located and is valid.
    return True, ("Matches alpha–amino acid motif: an amino group directly bound to an alpha carbon that connects "
                  "to a terminal carboxyl group, with no peptide bonds detected and acceptable core connectivity.")

    
# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "[H][C@@](\\N=C1CC(O)(CO)CC(NCC(O)=O)=C\\1OC)([C@H](C)O)C(O)=O",  # porphyra-334
        "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(O)=O",            # L-thyroxine
        "N[C@H](CCOc1ccc(cc1)C(=N\\O)\\C(=O)N[C@H]1CN([C@@H](C(O)=O)c2ccc(O)cc2)C1=O)C(O)=O",  # nocardicin A
        "O1C(C1C(=O)N)C(=O)NC[C@H](N)C(O)=O",                             # (2S)-3-(3-Carbamoyloxiranylcarbonylamino)-2-aminopropanoic acid
        "OC1CNC(C1)C(O)=O",                                               # 4-hydroxyproline
        "N[C@@H](C[C@@H]1CCC(=O)[C@@H]2O[C@H]12)C(O)=O",                  # anticapsin
        "N[C@@H](Cc1ccc(O)c(c1)[N+]([O-])=O)C(O)=O",                       # 3-nitro-L-tyrosine
        "NC(=O)CC[C@H](Nc1c(cc(cc1[N+]([O-])=O)[N+]([O-])=O)[N+]([O-])=O)C(O)=O",  # Tnp-Gln
        "CSCCCCCCC(N)C(O)=O",                                             # tetrahomomethionine
        "N[C@@H](CCSC[C@H]1OC(O)[C@H](O)[C@@H]1O)C(O)=O",                 # S-(5-deoxy-D-ribos-5-yl)-L-homocysteine
    ]
    for smi in test_smiles:
        result, reason = is_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")