"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: N-acyl-L-alpha-amino acid
Definition: Any L-alpha-amino acid carrying an N-acyl substituent.
An L-alpha-amino acid has a chiral alpha‐carbon bonded to an amino group and a carboxyl group.
The N-acyl substituent is defined by an acyl group (R–C(=O)–) attached directly to the amino nitrogen,
and this substituent must not be part of a peptide backbone.
"""

from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl L-alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is an N-acyl L-alpha-amino acid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1: Look for an L-alpha-amino acid backbone ---
    # We search for a chiral (C@ or C@@) carbon with an amino (N) group and a carboxyl (C(=O)[O,OH,-]) group.
    # Allowing for either neutral acid or deprotonated acid.
    aa_smarts = [
        "[C@H](N)(C(=O)[O;H,-])",
        "[C@@H](N)(C(=O)[O;H,-])"
    ]
    aa_matches = []
    for sm in aa_smarts:
        patt = Chem.MolFromSmarts(sm)
        aa_matches.extend(mol.GetSubstructMatches(patt))
    if not aa_matches:
        return False, "No L-alpha-amino acid backbone found"
    
    # Helper: Check if an atom (candidate acyl carbon) has an immediate carbonyl oxygen.
    def has_direct_carbonyl(carbon):
        # Look for a double bond (order==DOUBLE) to an oxygen
        for bond in carbon.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(carbon)
                if nbr.GetAtomicNum() == 8:
                    return True
        return False

    # Helper: Determine if the acyl group is part of a peptide bond.
    # In a typical peptide bond, the carbonyl carbon is linked not only to the amino nitrogen but 
    # also to another chiral carbon (alpha carbon from the previous residue).
    def is_peptide_bond(acyl_carbon, attached_nitrogen):
        for nbr in acyl_carbon.GetNeighbors():
            # Skip the attached nitrogen and any oxygen (from the carbonyl)
            if nbr.GetIdx() == attached_nitrogen.GetIdx():
                continue
            if nbr.GetAtomicNum() == 8:
                continue
            # If any other neighbor is a chiral carbon, suspect a peptide-bond connection.
            if nbr.GetAtomicNum() == 6 and nbr.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                return True
        return False

    # --- Step 2: For each candidate L-alpha-amino acid unit, check for an N-acyl substitution ---
    for match in aa_matches:
        # In our SMARTS, the chiral alpha carbon is the first atom of the match.
        alpha_idx = match[0]
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        # Find nitrogen atoms directly attached to the alpha carbon (the amino group).
        amino_nitrogens = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() == 7]
        if not amino_nitrogens:
            continue  # Shouldn't happen because our SMARTS required an N, but be safe.
        for amine in amino_nitrogens:
            # For each nitrogen, check its other substituents (neighbors except the alpha carbon).
            for substituent in amine.GetNeighbors():
                if substituent.GetIdx() == alpha_idx:
                    continue
                # We expect the acyl group to be a carbon (atomic number 6)
                if substituent.GetAtomicNum() != 6:
                    continue
                # Check that the bond between the amine and the substituent is a single bond.
                bond = mol.GetBondBetweenAtoms(amine.GetIdx(), substituent.GetIdx())
                if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                    continue
                # Now check that this carbon bears a carbonyl (has a double-bonded oxygen).
                if not has_direct_carbonyl(substituent):
                    continue
                # To avoid flagging peptide bonds, check if this acyl group appears connected to another chiral carbon.
                if is_peptide_bond(substituent, amine):
                    continue
                # Found a proper N-acyl group.
                return True, "Contains L-alpha-amino acid backbone with acylated amino group"
    
    return False, "Found L-alpha-amino acid backbone, but amino group is not acylated"

# When run as a script you can test with an example.
if __name__ == "__main__":
    # Example test: N-acetyl-L-aspartic acid
    example_smiles = "CC(=O)N[C@@H](CC(O)=O)C(O)=O"
    result, reason = is_N_acyl_L_alpha_amino_acid(example_smiles)
    print("Result:", result)
    print("Reason:", reason)