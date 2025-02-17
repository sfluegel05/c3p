"""
Classifies: CHEBI:83820 non-proteinogenic amino acid
"""
#!/usr/bin/env python
"""
Classifies: Non-proteinogenic amino acid 
Definition: Any amino-acid that is not naturally encoded in the genetic code of any organism.

The function is_non_proteinogenic_amino_acid takes a SMILES string as input and returns:
    (bool, str) -> (True, reason) if the molecule appears to be a free amino acid 
                     that is NOT one of the canonical proteinogenic ones,
                  (False, reason) otherwise.
                  
The algorithm applies several heuristics:
  1. Reject invalid SMILES.
  2. Explicitly add hydrogens so that the alpha-carbon (which must have at least one hydrogen)
     is correctly evaluated.
  3. Reject molecules that contain an internal peptide bond (SMARTS: "C(=O)N[C]").
  4. Identify a free carboxyl group (using the pattern "[CX3](=O)[OX1,OX2]") and require exactly one.
  5. Among the neighbors of the carboxyl carbon, look for a candidate “alpha‐carbon” that is:
       • sp3,
       • carries at least one hydrogen,
       • is attached to at least one nitrogen atom (the amino group),
       • and is not itself part of an amide bond.
  6. Finally, compute the InChIKey and if it matches one of the 20 canonical proteinogenic amino acids,
     then reject the molecule.
     
Note: This heuristic will not be perfect. Some non‐proteinogenic amino acids that lack a typical alpha‐H
(e.g. alpha‐methylated variants) might be missed, and some small peptides or amino acid derivatives might
falsely trigger the motif. This is one attempt to balance false positives versus false negatives.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

# Pre-compute InChIKeys for the 20 canonical proteinogenic amino acids.
_PROTEINOGENIC_SMILES = [
    "NCC(=O)O",                         # Glycine
    "CC(N)C(=O)O",                       # Alanine
    "CC(C)[C@H](N)C(=O)O",                # Valine
    "CCC[C@H](N)C(=O)O",                  # Leucine
    "CC[C@H](C)[C@H](N)C(=O)O",           # Isoleucine
    "C1CC[C@@H](NC1)C(=O)O",              # Proline
    "N[C@@H](Cc1ccccc1)C(=O)O",           # Phenylalanine
    "N[C@@H](Cc1ccc(O)cc1)C(=O)O",         # Tyrosine
    "N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O",    # Tryptophan
    "N[C@@H](CO)C(=O)O",                  # Serine
    "N[C@@H]([C@@H](O)C)C(=O)O",          # Threonine
    "N[C@@H](CS)C(=O)O",                  # Cysteine
    "N[C@@H](CCSC)C(=O)O",                # Methionine
    "N[C@@H](CC(=O)O)C(=O)O",             # Aspartic acid
    "N[C@@H](CCC(=O)O)C(=O)O",            # Glutamic acid
    "N[C@@H](CC(=O)N)C(=O)O",             # Asparagine
    "N[C@@H](CCC(=O)N)C(=O)O",            # Glutamine
    "N[C@@H](CCCCN)C(=O)O",               # Lysine
    "N[C@@H](CCCNC(=N)N)C(=O)O",          # Arginine
    "N[C@@H](Cc1cnc[nH]1)C(=O)O",          # Histidine
]
_PROTEINOGENIC_INCHIKEYS = set()
for smi in _PROTEINOGENIC_SMILES:
    mol = Chem.MolFromSmiles(smi)
    if mol is not None:
        ik = Chem.MolToInchiKey(mol)
        _PROTEINOGENIC_INCHIKEYS.add(ik)

def is_non_proteinogenic_amino_acid(smiles: str):
    """
    Determines whether the SMILES string represents a free amino acid that is
    non-proteinogenic (i.e. not one of the 20 canonical amino acids).
    
    The function checks:
      - That the molecule parses and explicit hydrogens are added.
      - That (a) there is exactly one free carboxyl group and 
           (b) at least one candidate alpha-carbon is found that is sp3, has at least one hydrogen,
               and is connected to at least one nitrogen (free amino group).
      - That no internal peptide bond (amide bond between amino acid units) is found.
      - That the molecule’s InChIKey is not in the set of canonical proteinogenic amino acids.
      
    Args:
        smiles (str): SMILES string.
        
    Returns:
        (bool, str): (True, reason) if classified as a non-proteinogenic amino acid;
                     (False, reason) otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens (needed to check hydrogen counts on the alpha-carbon)
    mol = Chem.AddHs(mol)
    
    # Reject if any peptide bond (internal amide linkage) is found.
    # This SMARTS looks for a C(=O)N[C] fragment. Even one such match is taken
    # as evidence of a peptide or dipeptide.
    peptide_bb = Chem.MolFromSmarts("C(=O)N[C]")
    if mol.HasSubstructMatch(peptide_bb):
        return False, "Molecule appears to be part of a peptide (has internal amide bonds)"
    
    # Identify free carboxyl group(s) using a SMARTS.
    # This pattern looks for a carbonyl where the carbon is in a carboxyl group.
    carboxyl_smarts = Chem.MolFromSmarts("[CX3](=O)[OX1,OX2]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_smarts)
    if len(carboxyl_matches) != 1:
        return False, f"Expected exactly one free carboxyl group, found {len(carboxyl_matches)}"
    
    # Assume the first matched atom is the carboxyl carbon.
    carboxyl_atom_idx = carboxyl_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_atom_idx)
    
    # Try to identify an alpha carbon candidate:
    # It should be directly bonded to the carboxyl carbon,
    # be an sp3 carbon with at least one explicit hydrogen,
    # and have at least one neighbor nitrogen that is not in an amide bond.
    alpha_candidate_found = False
    for nbr in carboxyl_atom.GetNeighbors():
        # We expect the alpha carbon to be carbon and not oxygen.
        if nbr.GetAtomicNum() != 6:
            continue
        # Check hybridization (sp3 preferred)
        if nbr.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        # Check that this carbon has at least one hydrogen (explicitly present)
        # Using GetTotalNumHs(includeNeighbors=True) after AddHs.
        if nbr.GetTotalNumHs() < 1:
            continue
        # Now, check that among its neighbors (other than the carboxyl carbon) 
        # there is at least one nitrogen that is not part of an amide bond.
        has_free_amino = False
        for subnbr in nbr.GetNeighbors():
            if subnbr.GetIdx() == carboxyl_atom.GetIdx():
                continue
            if subnbr.GetAtomicNum() == 7:
                # Check if the N is not part of an amide linkage.
                # If the nitrogen is bonded to a carbonyl which in turn is bonded to something else,
                # we treat it as being in a peptide bond.
                is_amide = False
                for nnbr in subnbr.GetNeighbors():
                    if nnbr.GetAtomicNum() == 6:
                        for o in nnbr.GetNeighbors():
                            if o.GetAtomicNum() == 8 and o.GetSymbol() == "O":
                                # Found a C=O near the N.
                                is_amide = True
                                break
                    if is_amide:
                        break
                if not is_amide:
                    has_free_amino = True
                    break
        if has_free_amino:
            alpha_candidate_found = True
            break
    if not alpha_candidate_found:
        return False, "No valid alpha-carbon with a free amino group found"

    # Compute the InChIKey of the (hydrogen-added) molecule.
    try:
        inchi_key = Chem.MolToInchiKey(mol)
    except Exception as e:
        return False, f"Failed to generate InChIKey: {e}"
    
    if inchi_key in _PROTEINOGENIC_INCHIKEYS:
        return False, "Matches a canonical proteinogenic amino acid"
    
    return True, "Contains a free carboxyl group and a candidate alpha-carbon with a free amino group, and does not match any canonical proteinogenic amino acid"
    
# Example usage:
# test_smiles = [
#    "[C@H]1(C(N(C1)[C@H](C=2C=CC(=CC2)O)C(=O)O)=O)NC([C@@H](C3=CC=C(C=C3)OCC[C@H](C(O)=O)N)=O", # nocardicin C (should be non-proteinogenic)
#    "NCCCC[C@H](N)CC(O)=O",  # (3S)-3,7-diaminoheptanoic acid (non-proteinogenic, target to be True)
#    "NCC(=O)O",             # Glycine (proteinogenic, should return False)
#    "FC(N)(C)C(O)=O",        # Aminofluoropropionic acid (false positive in previous attempt)
# ]
# for smi in test_smiles:
#     result, reason = is_non_proteinogenic_amino_acid(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")