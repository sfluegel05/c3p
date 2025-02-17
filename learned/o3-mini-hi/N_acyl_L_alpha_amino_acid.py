"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: N-acyl-L-alpha-amino acid
Definition: Any L-alpha-amino acid carrying an N-acyl substituent.
An L-alpha-amino acid has a chiral alpha‐carbon bonded to an amino group and a free carboxyl group.
The N-acyl substituent is defined as an acyl group (R–C(=O)–) attached directly to the nitrogen,
which must not be part of a peptide bond.
"""

from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl L-alpha-amino acid based on its SMILES string.

    The function first attempts to detect an L-alpha-amino acid unit by matching a chiral
    carbon attached to an amino group and a carboxyl group. To avoid peptide chains we
    require that the carboxyl group is “free” (i.e. the carboxyl carbon is not bonded to any nitrogen
    other than that from the backbone). Then, for the amino group bonded to that chiral carbon,
    we search for a substituent that meets the criteria for an N-acyl group: a single-bonded carbon
    that itself is doubly bound to an oxygen, and that is not connected (besides the carbonyl oxygen)
    to another chiral center typical of a peptide bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a valid N-acyl L-alpha-amino acid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS to capture L-alpha-amino acid backbone:
    # It looks for a chiral carbon ([C@H] or [C@@H]) connected to an amino group (N)
    # and a carboxyl group (C(=O)[O;H,-]), which should be free (not connected to another N).
    aa_smarts = [
        "[C@H](N)(C(=O)[O;H,-])",
        "[C@@H](N)(C(=O)[O;H,-])"
    ]
    backbone_matches = []
    for sm in aa_smarts:
        patt = Chem.MolFromSmarts(sm)
        if patt is None:
            continue
        backbone_matches.extend(mol.GetSubstructMatches(patt))
        
    if not backbone_matches:
        return False, "No L-alpha-amino acid backbone found"

    # Helper function to check if a carbon atom is an acyl carbon
    # i.e. it has a carbonyl oxygen (double-bonded O).
    def has_carbonyl(carbon):
        for bond in carbon.GetBonds():
            # check for double bond to oxygen
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(carbon)
                if nbr.GetAtomicNum() == 8:
                    return True
        return False

    # Helper: determine if the acyl carbon appears to be part of a peptide bond.
    # In a peptide bond, the acyl carbon is connected to a chiral (alpha) carbon from a neighboring residue.
    def seems_peptide_bond(acyl_carbon, attached_N):
        # Check all heavy neighbors except the one we came from (attached_N) and carbonyl oxygens.
        for nbr in acyl_carbon.GetNeighbors():
            if nbr.GetIdx() == attached_N.GetIdx():
                continue
            # Skip oxygen from carbonyl (double bond)
            if nbr.GetAtomicNum() == 8:
                continue
            # If any neighbor is a chiral carbon with specified chirality, suspect peptide bond.
            if nbr.GetAtomicNum() == 6 and nbr.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                return True
        return False

    # When dealing with peptides, the carboxyl group is amide-bound.
    # We require the carboxyl carbon of the backbone to have NO nitrogen neighbor.
    def is_free_carboxyl(carboxyl_carbon):
        for nbr in carboxyl_carbon.GetNeighbors():
            # Skip oxygen neighbors (the carbonyl and hydroxyl oxygen)
            if nbr.GetAtomicNum() == 8:
                continue
            if nbr.GetAtomicNum() == 7:
                return False
        return True

    # For each identified amino acid backbone match, check:
    for match in backbone_matches:
        # Our SMARTS was defined as [C](N)(C(=O)[O;H,-]), so we expect three atoms;
        # assume the first atom in the match is the chiral alpha carbon.
        alpha_idx = match[0]
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        # In the match, we expect a nitrogen and a carboxyl-carbon.
        # Distinguish the two by their atomic number and bond environment.
        neighbor_atoms = alpha_atom.GetNeighbors()
        candidate_N = None
        candidate_C_carboxyl = None
        for nbr in neighbor_atoms:
            if nbr.GetAtomicNum() == 7:
                candidate_N = nbr
            elif nbr.GetAtomicNum() == 6:
                # Check if this carbon atom has a carbonyl pattern (i.e., is part of COOH)
                if has_carbonyl(nbr):
                    candidate_C_carboxyl = nbr
        if candidate_N is None or candidate_C_carboxyl is None:
            continue  # not matching a complete AA unit

        # Check that the carboxyl group is free (not involved in peptide bond)
        if not is_free_carboxyl(candidate_C_carboxyl):
            continue  # likely a di-/oligo-peptide

        # Now, for the amine candidate (the N attached to alpha),
        # look at its substituents other than the alpha carbon.
        found_acyl = False
        for nbr in candidate_N.GetNeighbors():
            if nbr.GetIdx() == alpha_atom.GetIdx():
                continue
            # We are looking for a carbon that is connected via a single bond:
            if nbr.GetAtomicNum() != 6:
                continue
            bond = mol.GetBondBetweenAtoms(candidate_N.GetIdx(), nbr.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            acyl_carbon = nbr
            # Check that this carbon has a double-bonded oxygen (i.e. is carbonyl carbon)
            if not has_carbonyl(acyl_carbon):
                continue
            # To exclude peptide bonds (where the acyl carbon is part of the backbone), 
            # check that it is not additionally connected to a chiral carbon.
            if seems_peptide_bond(acyl_carbon, candidate_N):
                continue
            # It appears we found an N-acyl substituent.
            found_acyl = True
            break

        if found_acyl:
            return True, "Contains L-alpha-amino acid backbone with acylated amino group"

    return False, "Found L-alpha-amino acid backbone with free carboxyl group, but amino group is not acylated"

# When run as a script you can test with an example.
if __name__ == "__main__":
    # Example test: N-acetyl-L-aspartic acid
    example_smiles = "CC(=O)N[C@@H](CC(O)=O)C(O)=O"
    result, reason = is_N_acyl_L_alpha_amino_acid(example_smiles)
    print("Result:", result)
    print("Reason:", reason)