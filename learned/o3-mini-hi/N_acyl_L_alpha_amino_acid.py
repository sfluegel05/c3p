"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: N-acyl-L-alpha-amino acid.
Definition: Any L-alpha-amino acid carrying an N-acyl substituent.
An L-alpha-amino acid is defined as having a chiral α–carbon bonded to a primary amine and a free carboxyl group.
The N-acyl substituent is an acyl group (R–C(=O)–) that is attached directly to a nitrogen (which may be the backbone’s amino group or a side‐chain nitrogen),
with the acyl carbon not being part of an inter-residue peptide bond.
This implementation first identifies a unique “free” amino acid backbone and then looks for a substituent that shows the N–C(=O) pattern.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines whether a molecule (provided as a SMILES string) qualifies as an N-acyl L-α-amino acid.

    The function first attempts to locate a unique L-α-amino acid backbone. We require:
      • A chiral α–carbon (C@ or C@@) bonded to a primary amine (–NH2) and
      • a free carboxyl group (–C(=O)O or –C(=O)[O-]) that is not engaged in a peptide bond.
    Then, the molecule is scanned for any nitrogen that bears an acyl substituent. An acyl substituent is defined
    as a carbon bound via a single bond to nitrogen that carries at least one double-bonded oxygen (i.e. a carbonyl).
    To avoid mis‐identifying an amide bond within a peptide chain as an “acyl substituent,” we skip the N–C connection 
    if that carbon is already the carboxyl carbon of the backbone.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an N-acyl L-α-amino acid, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # ---- Step 1. Locate exactly one free L-α-amino acid backbone.
    # We use a stricter SMARTS that requires:
    #   - a chiral carbon ([C@H] or [C@@H])
    #   - bonded to a nitrogen (the amino group)
    #   - bonded to a carboxyl group carbon: C(=O)[O;H,-]
    # Also, we will later check that the carboxyl carbon does not have any nitrogen neighbors
    # (which would indicate it is involved in a peptide bond).
    backbone_smarts_list = [
        "[C@H](N)(C(=O)[O;H1,-1])",
        "[C@@H](N)(C(=O)[O;H1,-1])",
        "[C@H](N)(C(=O)O)",
        "[C@@H](N)(C(=O)O)"
    ]
    backbone_matches = []
    match_detail = None  # record details (indices) of accepted backbone
    for smarts in backbone_smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        matches = mol.GetSubstructMatches(patt)
        for m in matches:
            backbone_matches.append(m)
    # Remove duplicates by sorting tuples
    backbone_matches = list(set(tuple(sorted(match)) for match in backbone_matches))
    
    if not backbone_matches:
        return False, "No L-alpha-amino acid backbone found"
    if len(backbone_matches) > 1:
        # This likely indicates a peptide chain rather than a single amino acid.
        return False, "Multiple amino acid backbones detected (likely a peptide chain)"
    
    backbone_atoms = backbone_matches[0]
    # In our SMARTS, the matching atom order is not guaranteed.
    # We now try to determine which atom is:
    #   alpha carbon: must be carbon, chiral and bonded to at least one nitrogen and one carbon.
    #   backbone amino nitrogen: atomic num 7.
    #   carboxyl carbon: a carbon that has a C=O bond, and must not be bound to any additional nitrogen.
    alpha_idx = None
    amine_idx = None
    carboxyl_idx = None
    for idx in backbone_atoms:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() == 6 and atom.GetChiralTag() != rdchem.ChiralType.CHI_UNSPECIFIED:
            # candidate for alpha carbon
            alpha_idx = idx
        elif atom.GetAtomicNum() == 7:
            # candidate for amine nitrogen.
            amine_idx = idx
        # Try to identify carboxyl carbon: a carbon with at least one double-bonded oxygen neighbor.
        if atom.GetAtomicNum() == 6:
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bond and bond.GetBondType() == rdchem.BondType.DOUBLE:
                        carboxyl_idx = idx
                        break
    if alpha_idx is None or amine_idx is None or carboxyl_idx is None:
        return False, "Incomplete amino acid backbone identified"

    # Check the carboxyl carbon to be "free" (i.e. not forming a peptide bond).
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 7:
            # If the carboxyl carbon has any nitrogen neighbor (beyond its carbonyl O),
            # it may be involved in a peptide bond.
            return False, "Carboxyl group is not free (likely involved in a peptide bond)"
    
    # ---- Step 2. Search for an N-acyl substituent.
    # We look for any nitrogen atom that is bound (by a single bond) to a carbon that has a double-bonded oxygen.
    acyl_found = False
    acyl_atom_info = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        # For the backbone amine nitrogen, ignore the bond to the backbone alpha carbon.
        for nbr in atom.GetNeighbors():
            # Only consider carbon neighbors.
            if nbr.GetAtomicNum() != 6:
                continue
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if not bond or bond.GetBondType() != rdchem.BondType.SINGLE:
                continue
            # If this carbon is the carboxyl carbon from our backbone, skip it.
            if nbr.GetIdx() == carboxyl_idx:
                continue
            # Check if this carbon (candidate acyl carbon) has at least one double-bonded oxygen.
            has_carbonyl = False
            for sub in nbr.GetNeighbors():
                if sub.GetAtomicNum() == 8:
                    bond_sub = mol.GetBondBetweenAtoms(nbr.GetIdx(), sub.GetIdx())
                    if bond_sub and bond_sub.GetBondType() == rdchem.BondType.DOUBLE:
                        has_carbonyl = True
                        break
            if has_carbonyl:
                # We have found an N–C(=O) connection that appears as an acyl substituent.
                acyl_found = True
                acyl_atom_info = (atom.GetIdx(), nbr.GetIdx())
                break
        if acyl_found:
            break

    if acyl_found:
        return True, "Contains L-alpha-amino acid backbone with acylated amino group"
    else:
        return False, "Found L-alpha-amino acid backbone, but no N-acyl substituent detected"

# When invoked as a script, run an example test.
if __name__ == "__main__":
    # Example test: one of the provided N-acetyl-L-aspartic acid examples.
    example_smiles = "CC(=O)N[C@@H](CC(O)=O)C(O)=O"
    result, reason = is_N_acyl_L_alpha_amino_acid(example_smiles)
    print("Result:", result)
    print("Reason:", reason)