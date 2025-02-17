"""
Classifies: CHEBI:33704 alpha-amino acid
"""
#!/usr/bin/env python
"""
Classifies: Alpha–amino acid
Definition: An amino acid in which the amino group is located on the carbon atom 
at the position alpha to the carboxy group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid(smiles: str):
    """
    Determines whether a molecule (given by its SMILES) is a stand-alone alpha–amino acid.
    
    The algorithm:
      1. Parses the SMILES and adds explicit hydrogens.
      2. Rejects molecules that contain a peptide (amide) bond using a SMARTS.
         (This prevents di-/oligo–peptides from matching.)
      3. Rejects molecules that are unusually large (heuristic heavy atom count threshold).
      4. Identifies a free carboxyl group using a SMARTS pattern:
         The pattern requires a carboxyl carbon (sp2, with a double bond to one oxygen and a single bond 
         to a hydroxyl (or deprotonated –O–) oxygen). We then require that this carboxyl carbon has exactly 
         one neighbor that is not oxygen (the candidate alpha–carbon).
      5. For the candidate alpha–carbon, verify that it is sp3 (or nearly tetrahedral) and that it is directly 
         bonded to at least one nitrogen. For non–ring nitrogen we require at least two attached hydrogens.
      6. Require that exactly one motif is found.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): A tuple containing the classification (True if an alpha–amino acid) and an explanation.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for peptide (amide) bonds. (A peptide bond connects two amino acids)
    # Here we look for "C(=O)N[C]" i.e. an amide carbonyl with a nitrogen that is bound further.
    peptide_smarts = Chem.MolFromSmarts("C(=O)N[C]")
    if mol.HasSubstructMatch(peptide_smarts):
        return False, "Contains peptide bond(s) indicating a peptide rather than a single amino acid"
    
    # Reject molecules that are too large to be standard standalone amino acids.
    if mol.GetNumHeavyAtoms() > 50:
        return False, "Molecule size too large to be a simple amino acid"
    
    # Add explicit hydrogens to reliably check for hydrogen counts on nitrogen.
    molH = Chem.AddHs(mol)
    
    # SMARTS for a free carboxyl group.
    # This matches a carboxyl carbon: it is sp2, double-bonded to an oxygen and single-bonded to an OH (or O-)
    # Note: We do not force the charge on the oxygen so that both protonated and deprotonated forms are caught.
    carboxyl_smarts = Chem.MolFromSmarts("[CX3](=O)[O;H1,-]")
    carboxyl_matches = molH.GetSubstructMatches(carboxyl_smarts)
    
    motif_count = 0
    reasons = []
    # For every carboxyl match:
    # The SMARTS returns a tuple of indices (c-o-hydroxy, carbonyl oxygen, hydroxyl oxygen).
    # We want to get the carboxyl carbon and then examine its non-oxygen neighbor.
    for match in carboxyl_matches:
        carboxyl_c = molH.GetAtomWithIdx(match[0])
        # Get neighbors of carboxyl carbon that are not oxygen.
        alpha_candidates = [nbr for nbr in carboxyl_c.GetNeighbors() if nbr.GetAtomicNum() != 8]
        if len(alpha_candidates) != 1:
            # To be a free terminal carboxyl group in an amino acid, the carboxyl carbon should have 
            # exactly one non-oxygen neighbor (the alpha carbon).
            continue
        alpha_c = alpha_candidates[0]
        # (Optional) Check that the carboxyl carbon is degree 3 (two oxygens and one alpha candidate).
        if len(carboxyl_c.GetNeighbors()) != 3:
            continue
        
        # Now, for the candidate alpha carbon, verify it is a carbon and likely tetrahedral (sp3).
        if alpha_c.GetAtomicNum() != 6:
            continue
        # In a typical amino acid, the alpha carbon is sp3 (or at least not double-bonded).
        # We allow ring atoms (e.g. in proline derivatives) but otherwise expect no double bonds.
        if alpha_c.GetHybridization().name not in ["SP3", "UNSPECIFIED"]:
            continue
        
        # Look for an amino group attached to the alpha carbon.
        nitrogen_found = False
        for neighbor in alpha_c.GetNeighbors():
            if neighbor.GetAtomicNum() == 7:
                # If the nitrogen is not in a ring (non–proline-like), check that it carries at least two hydrogens.
                # In many zwitterionic forms the H count may be adjusted, but we at least try to filter out quaternary nitrogens.
                if not neighbor.IsInRing():
                    if neighbor.GetTotalNumHs() < 2:
                        continue
                nitrogen_found = True
                break
        if not nitrogen_found:
            continue
        
        # If we get here, we found one valid motif.
        motif_count += 1

    if motif_count == 1:
        return True, ("Matches alpha–amino acid motif: a free carboxyl group attached to an alpha carbon that "
                      "bears an amino group (with appropriate hydrogen count) and no peptide bonds are present")
    elif motif_count > 1:
        return False, f"Multiple potential alpha–amino acid motifs found ({motif_count}); not a simple amino acid"
    else:
        return False, "No alpha–amino acid motif found"


# Example usage:
# You can test on one or more SMILES from the provided list:
if __name__ == "__main__":
    test_smiles = [
       "[H][C@@](\\N=C1CC(O)(CO)CC(NCC(O)=O)=C\\1OC)([C@H](C)O)C(O)=O",  # porphyra-334, known true
       "N[C@@H](Cc1cc(I)c(Oc2cc(I)c(O)c(I)c2)c(I)c1)C(O)=O",            # L-thyroxine, known true
       "NC1(CC1)C(O)=O",                                                # 1-aminocyclopropanecarboxylic acid, expected true
       "[O-]C(=O)[C@H]1CNCC[NH2+]1",                                     # (R)-piperazine-2-carboxylic acid zwitterion, known false positive
    ]
    for smi in test_smiles:
        result, reason = is_alpha_amino_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")