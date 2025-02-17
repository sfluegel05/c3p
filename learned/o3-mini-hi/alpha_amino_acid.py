"""
Classifies: CHEBI:33704 alpha-amino acid
"""
#!/usr/bin/env python
"""
Classifies: Alpha-amino acid
Definition: An amino acid in which the amino group is located on the carbon atom at the position alpha to the carboxy group.
Revised approach:
  • Reject molecules that contain peptide (amide) bonds.
  • Reject molecules that are unusually large (likely peptides or complex conjugates).
  • Identify a free carboxyl group by finding a carbon that is double‐bonded to one oxygen and single‐bonded to an –OH (or deprotonated –O–).
    Its third neighbor (excluding the oxygens) is considered the candidate “alpha carbon.”
  • Check that the candidate “alpha carbon” is attached (by a single bond) to at least one nitrogen.
    To reduce false positives, if the nitrogen is not in a ring then it must display at least two hydrogens.
  • Finally, require that exactly one such motif is found.
Note: This heuristic is still not perfect, but it reduces false positives from larger and multi–functional molecules.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is a stand-alone alpha–amino acid based on its SMILES string.
    The algorithm:
      1. Rejects molecules containing peptide (amide) bonds.
      2. Adds explicit hydrogens and rejects very large molecules.
      3. Searches for a free carboxyl group defined as a carboxyl carbon that is
         double–bonded to one oxygen and single–bonded to one oxygen (–OH or –O–).
         The third bond from that carboxyl carbon should be to an “alpha–carbon.”
      4. The candidate alpha–carbon must have at least one bonded nitrogen.
         For non–ring (i.e. non–proline) nitrogens, at least two hydrogens are required.
      5. To be accepted as an isolated alpha–amino acid, exactly one such motif must be present.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as (stand-alone) alpha–amino acid, False otherwise.
        str: Explanation of the decision.
    """
    # Parse SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject if any peptide (amide) bonds are present.
    # (This helps avoid di-/oligopeptides even if the amino terminus is free.)
    amide_bond_smarts = Chem.MolFromSmarts("C(=O)N")
    if mol.HasSubstructMatch(amide_bond_smarts):
        return False, "Contains peptide bond(s) indicating a peptide rather than a single amino acid"
    
    # Reject molecules that are too large to be typical standalone amino acids.
    # (Threshold chosen heuristically; most natural amino acids are small.)
    if mol.GetNumHeavyAtoms() > 50:
        return False, "Molecule size too large to be a single amino acid"
    
    # Add explicit hydrogens so we can verify hydrogen counts on amino nitrogens.
    molH = Chem.AddHs(mol)
    
    motif_count = 0
    # For each atom, try to identify a free carboxyl group.
    # A free carboxyl carbon is defined as a carbon atom that:
    #   - Has a double bond to one oxygen (carbonyl oxygen).
    #   - Has a single bond to one oxygen that is part of an –OH (or deprotonated –O–, but may show H attached).
    #   - Has exactly one other neighbor (the candidate alpha–carbon).
    for atom in molH.GetAtoms():
        # Look only at carbon atoms.
        if atom.GetAtomicNum() != 6:
            continue
        neighbors = atom.GetNeighbors()
        # We expect exactly three neighbors for an ideal carboxyl carbon.
        if len(neighbors) != 3:
            continue
        carbonyl_oxygen = None
        hydroxyl_oxygen = None
        alpha_candidate = None
        # Examine neighbors.
        for nbr in neighbors:
            bond = molH.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # If neighbor is oxygen
            if nbr.GetAtomicNum() == 8:
                # Check bond order (using RDKit bond order char)
                bond_order = bond.GetBondTypeAsDouble()
                if bond_order == 2 and carbonyl_oxygen is None:
                    carbonyl_oxygen = nbr
                elif bond_order == 1 and hydroxyl_oxygen is None:
                    # For OH we require at least one hydrogen attached,
                    # though in deprotonated forms the explicit hydrogen might be missing.
                    # So we check either explicit H count > 0 or a negative charge.
                    if nbr.GetTotalNumHs() >= 1 or nbr.GetFormalCharge() < 0:
                        hydroxyl_oxygen = nbr
            else:
                # That neighbor is the candidate alpha–carbon.
                if alpha_candidate is None:
                    alpha_candidate = nbr
        # Must have both oxygens and an alpha candidate.
        if carbonyl_oxygen is None or hydroxyl_oxygen is None or alpha_candidate is None:
            continue
        # Now, check the alpha candidate:
        # It should be a carbon (for glycine or any amino acid; note that in proline the alpha carbon may be cyclic, which is allowed).
        if alpha_candidate.GetAtomicNum() != 6:
            continue
        # Look among alpha_candidate's neighbors for an amino nitrogen.
        nitrogen_found = False
        for nb in alpha_candidate.GetNeighbors():
            # We require that the bond is a single bond (should be the case for an amino group).
            if nb.GetAtomicNum() == 7:
                # Additional check: if the nitrogen is not in a ring (i.e. not proline)
                # then it should have at least 2 attached hydrogens.
                if not nb.IsInRing():
                    if nb.GetTotalNumHs() < 2:
                        continue
                nitrogen_found = True
                break
        if not nitrogen_found:
            continue
        
        # If we reached here, we found one instance of the motif.
        motif_count += 1

    if motif_count == 1:
        return True, ("Matches alpha–amino acid motif: a free carboxyl group attached to an alpha carbon that " +
                      "bears an amino group (with appropriate hydrogen count) and no peptide bonds present")
    elif motif_count > 1:
        return False, f"Multiple potential alpha–amino acid motifs found ({motif_count}); not a simple amino acid"
    else:
        return False, "No alpha–amino acid motif found"

# Example test cases (you may uncomment these to try the function):
# test_smiles = [
#     "[H][C@@](\\N=C1CC(O)(CO)CC(NCC(O)=O)=C\\1OC)([C@H](C)O)C(O)=O",  # porphyra-334, expected True
#     "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(O)=O",            # L-thyroxine, expected True
#     "NC1(CC1)C(O)=O",                                               # 1-aminocyclopropanecarboxylic acid, expected True
#     "OC1CNC(C1)C(O)=O",                                              # 4-hydroxyproline, expected True
#     "C(C([O-])=O)(CO)(CO)[NH3+]",                                    # alpha-(hydroxymethyl)serine zwitterion, false positive in previous run
#     "O=C(NCC(=O)N[C@@H](CO)C(O)=O)[C@@H](N)CC=1NC=NC1"               # A di/peptide, expected False
# ]
#
# for smi in test_smiles:
#     result, reason = is_alpha_amino_acid(smi)
#     print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")