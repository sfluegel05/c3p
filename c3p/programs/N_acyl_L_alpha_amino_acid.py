"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: N-acyl-L-alpha-amino acid
Definition: Any L-alpha-amino acid carrying an N-acyl substituent.
An L-alpha-amino acid is defined as having a chiral α–carbon bonded to an amino group and a free carboxyl group.
The N-acyl substituent is an acyl group (R–C(=O)–) that is attached directly to a nitrogen (which may be the backbone’s amino group or a side‐chain nitrogen),
with the acyl carbon not further attached to another chiral center (i.e. not forming a peptide bond).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is an N-acyl L-α-amino acid.
    
    First, this function looks for a unique L-α-amino acid backbone. We use several SMARTS to capture both protonated 
    and deprotonated forms of the carboxyl group. To avoid classifying peptides as single amino acids we require that 
    only one such backbone is present and that the carboxyl carbon is free (i.e. not bound to any nitrogen aside from those in the unit itself).
    
    Then the function scans the molecule for any nitrogen that carries an acyl substituent. An acyl substituent here is 
    defined as a carbon (bound by a single bond to the nitrogen) that itself is carbonyl‐like (has at least one double‐bonded oxygen)
    and is not part of a peptide‐like linkage.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an N-acyl L-α-amino acid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for an L-α-amino acid backbone.
    # These patterns capture a chiral carbon bonded to a primary amine (N) and a carboxyl group (as COOH or COO-).
    backbone_smarts_list = [
        "[C@H](N)(C(=O)[O;H,-])",
        "[C@@H](N)(C(=O)[O;H,-])",
        "[C@H](N)(C(=O)O)",
        "[C@@H](N)(C(=O)O)"
    ]
    backbone_matches = []
    for smarts in backbone_smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            backbone_matches.extend(mol.GetSubstructMatches(patt))
    # Remove any duplicate matches (order of atom indices may differ)
    backbone_matches = list(set(backbone_matches))
    
    if not backbone_matches:
        return False, "No L-alpha-amino acid backbone found"
    if len(backbone_matches) > 1:
        return False, "Multiple amino acid backbones detected (likely a peptide chain)"
    
    # For the unique backbone match, try to identify: a chiral (alpha) carbon, a nitrogen, and the carboxyl carbon.
    match = backbone_matches[0]
    alpha_idx = None
    n_idx = None
    carboxyl_idx = None
    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        # Look for the chiral alpha carbon (should have explicit chirality)
        if atom.GetAtomicNum() == 6 and atom.GetChiralTag() != rdchem.ChiralType.CHI_UNSPECIFIED:
            alpha_idx = idx
        elif atom.GetAtomicNum() == 7:
            # It could be a backbone amino nitrogen
            n_idx = idx
        elif atom.GetAtomicNum() == 6:
            # A candidate for the carboxyl carbon: look for a carbonyl pattern (attached by a double bond to oxygen)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                        carboxyl_idx = idx
                        break

    if alpha_idx is None or n_idx is None or carboxyl_idx is None:
        return False, "Incomplete amino acid backbone identified"
    
    # Check that the carboxyl carbon is "free": it should not be bound to any nitrogen (aside from those in its COOH group).
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 7:
            return False, "Carboxyl group is not free (likely involved in a peptide bond)"
    
    # Now search for any N-acyl substituent anywhere in the molecule.
    # An acyl substituent must be directly attached to a nitrogen via a single bond and consist of a carbon
    # that has at least one double-bonded oxygen. We also try to exclude cases where the acyl group is part of a peptide bond.
    acyl_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # nitrogen
            # Check all neighbors (other than those that are part of the backbone if needed)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() != 6:
                    continue
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
                    continue
                # Now check if 'nbr' (a candidate acyl carbon) has a double-bonded oxygen
                has_carbonyl = False
                for sub in nbr.GetNeighbors():
                    if sub.GetAtomicNum() == 8:
                        b_tmp = mol.GetBondBetweenAtoms(nbr.GetIdx(), sub.GetIdx())
                        if b_tmp and b_tmp.GetBondType() == rdchem.BondType.DOUBLE:
                            has_carbonyl = True
                            break
                if not has_carbonyl:
                    continue
                # To avoid peptide bonds we check that the candidate acyl carbon is not further attached 
                # to any chiral carbon (which would represent another backbone) aside from our current nitrogen.
                peptide_like = False
                for sub in nbr.GetNeighbors():
                    if sub.GetIdx() == atom.GetIdx():
                        continue
                    if sub.GetAtomicNum() == 6 and sub.GetChiralTag() != rdchem.ChiralType.CHI_UNSPECIFIED:
                        peptide_like = True
                        break
                if peptide_like:
                    continue
                # We found an acyl substituent
                acyl_found = True
                break
            if acyl_found:
                break

    if acyl_found:
        return True, "Contains L-alpha-amino acid backbone with acylated amino group"
    else:
        return False, "Found L-alpha-amino acid backbone, but no N-acyl substituent detected"

# When run as a script you can test with an example.
if __name__ == "__main__":
    # Example test: N-acetyl-L-aspartic acid (one of the examples provided)
    example_smiles = "CC(=O)N[C@@H](CC(O)=O)C(O)=O"
    result, reason = is_N_acyl_L_alpha_amino_acid(example_smiles)
    print("Result:", result)
    print("Reason:", reason)