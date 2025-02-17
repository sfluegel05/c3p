"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: N-acyl-L-alpha-amino acid 
Definition: Any L-alpha-amino acid carrying an N-acyl substituent.
An L-alpha-amino acid has a chiral alpha‐carbon bonded to an amino group and a carboxyl group.
The N-acyl substituent is defined by an acyl group (R–C(=O)–) attached directly to a nitrogen 
(which should not be part of a peptide backbone).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl L-alpha-amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is an N-acyl L-alpha-amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1: Look for an L-alpha-amino acid backbone ---
    # We require a chiral center (C@ or C@@) with an amino (N) substituent and a carboxyl group.
    # This pattern is a first approximation of [C@H](N)(C(=O)[O,-]) where the acid may be neutral or deprotonated.
    aa_smarts1 = Chem.MolFromSmarts("[C@H](N)(C(=O)[O;H,-])")
    aa_smarts2 = Chem.MolFromSmarts("[C@@H](N)(C(=O)[O;H,-])")
    matches1 = mol.GetSubstructMatches(aa_smarts1)
    matches2 = mol.GetSubstructMatches(aa_smarts2)
    aa_matches = matches1 + matches2
    if not aa_matches:
        return False, "No L-alpha-amino acid backbone found"
    
    # --- Helper: Check if a carbon atom exhibits a carbonyl group (i.e. double bond to oxygen) ---
    def has_carbonyl(carbon):
        # Return True if the carbon has at least one double bond to an oxygen atom.
        for bond in carbon.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(carbon)
                if nbr.GetAtomicNum() == 8:
                    return True
        return False

    # --- Helper: Determine whether the carbonyl carbon is part of a peptide bond.
    def appears_in_peptide_bond(acyl_carbon, acylated_nitrogen):
        """
        In a peptide bond the acyl carbon is attached not only to acylated_nitrogen but also to
        another chiral carbon (the alpha carbon of the preceding residue). If any neighboring atom
        (other than the acylated nitrogen and the double-bonded oxygen(s)) is chiral, we assume this
        carbonyl carbon is part of a peptide bond.
        """
        for nbr in acyl_carbon.GetNeighbors():
            # Exclude the amine already participating and the oxygen(s) in the carbonyl
            if nbr.GetIdx() == acylated_nitrogen.GetIdx():
                continue
            if nbr.GetAtomicNum() == 8:
                continue
            # Check chiral tag – if it is set (and not unspecified) we assume peptide connectivity.
            if nbr.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                return True
        return False

    # --- Step 2: For each found L-alpha-amino acid backbone, check for N-acyl substitution ---
    # We inspect nitrogen atoms directly attached to the alpha carbon.
    for match in aa_matches:
        # In the matched pattern, the first atom is the chiral alpha carbon.
        alpha_idx = match[0]
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        # Find all neighbor atoms that are nitrogen (atomic number 7)
        amino_nitrogens = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() == 7]
        if not amino_nitrogens:
            # Should not occur because our SMARTS required an N, but we check for safety.
            continue
        # For each amino nitrogen, inspect substituents (neighbors excluding the alpha carbon)
        for amine in amino_nitrogens:
            for sub in amine.GetNeighbors():
                if sub.GetIdx() == alpha_idx:
                    continue
                # Look for a carbon attached to the nitrogen
                if sub.GetAtomicNum() == 6:
                    # Check that this carbon has a double bond to oxygen (i.e. is a carbonyl carbon)
                    if has_carbonyl(sub):
                        # If the acyl carbon appears to be part of a peptide bond, skip it.
                        if appears_in_peptide_bond(sub, amine):
                            continue
                        # Otherwise we found an acyl group on the amino nitrogen.
                        return True, "Contains L-alpha-amino acid backbone with acylated amino group"
    # No suitable N-acyl substitution was found.
    return False, "Found L-alpha-amino acid backbone, but amino group is not acylated"

# When run as a script you can test with an example.
if __name__ == "__main__":
    # Example test: N-acetyl-L-aspartic acid
    example_smiles = "CC(=O)N[C@@H](CC(O)=O)C(O)=O"
    result, reason = is_N_acyl_L_alpha_amino_acid(example_smiles)
    print("Result:", result)
    print("Reason:", reason)