"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: Polyamine
Definition: Any organic amino compound that contains two or more amino groups.
We rule out peptide‐like compounds (which contain multiple amide bonds and typically only one “free” amino group)
and count as amino group any nitrogen that is not part of an amide.
Note: This is a heuristic approach using RDKit.
"""

from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.
    In addition to simply counting N atoms that are not in amide bonds,
    this function first checks whether the molecule appears to be a peptide (i.e. contains multiple amide bonds)
    and then counts any non-amide nitrogen atoms (including tertiary amines that may lack explicit H's).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a polyamine, False otherwise.
        str: A reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Check for peptide-like features.
    # We use a SMARTS pattern for an amide: a nitrogen directly bonded to a carbonyl carbon.
    amide_smarts = Chem.MolFromSmarts("[NX3][C](=[O])")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) >= 2:
        # If there are at least 2 amide bonds AND the molecule is not extremely small,
        # then it is likely a peptide (we want to prevent false positives).
        if mol.GetNumHeavyAtoms() > 10:
            return False, f"Molecule appears to be peptide-like ({len(amide_matches)} amide bonds)"
    
    # Step 2. Count potential amino groups.
    # We count a nitrogen as an amino group if:
    #  a) It is not directly involved in an amide bond (i.e. bonded to a carbon that has a double bond to oxygen)
    #  b) It is not aromatic within a heterocycle unless it is exocyclic (attached externally to an aromatic ring)
    #  c) We allow tertiary (and secondary, primary) aliphatic amines even if no explicit hydrogen is recorded.
    amino_count = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # not nitrogen
        
        # Check if the nitrogen is bonded to a carbon that is carbonyl (i.e. double-bonded to oxygen)
        is_amide = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # carbon neighbor
                for bond in neighbor.GetBonds():
                    # skip the bond back to the nitrogen itself
                    if bond.GetOtherAtom(neighbor).GetIdx() == atom.GetIdx():
                        continue
                    # if the neighbor carbon is double-bonded to oxygen, mark as amide
                    if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(neighbor).GetAtomicNum() == 8:
                        is_amide = True
                        break
            if is_amide:
                break
        if is_amide:
            continue  # do not count nitrogens that are directly in an amide environment
        
        # We now check ring membership and hybridization.
        # Do not count nitrogen that is aromatic and embedded in a heterocycle
        # (unless it is exocyclic; here we check if all its bonds are within a ring)
        if atom.GetIsAromatic():
            # if all bonds are in rings then we skip it (likely part of an aromatic heterocycle)
            if all(bond.IsInRing() for bond in atom.GetBonds()):
                continue
        
        # Exclude quaternary ammonium centers. (We allow tertiary amines even with 0 H)
        if atom.GetDegree() >= 4:
            continue
        
        # Otherwise, count this nitrogen.
        amino_count += 1

    # Decision: need at least two amino groups.
    if amino_count >= 2:
        return True, f"Contains {amino_count} amino groups."
    else:
        return False, f"Found only {amino_count} amino group(s), need at least 2 to be a polyamine."

# Example usage:
if __name__ == "__main__":
    # a few representative examples from the list
    test_smiles = [
        "CCNc1nc(N)nc(O)n1",                      # 4-amino-6-(ethylamino)-1,3,5-triazin-2-ol; expected polyamine
        "NCCCNCCNCCCN",                           # 3,2,3-tetramine; expected polyamine
        "NCCN",                                   # ethylenediamine; expected polyamine
        "CC(=O)NCCCCN",                           # N-acetylputrescine; expected polyamine
        "C1CN2CCN1CC2",                           # triethylenediamine; expected polyamine
        "O=C(N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O)...",  # Peptide-like (tripeptide) should be rejected
    ]
    for s in test_smiles:
        result, reason = is_polyamine(s)
        print(f"SMILES: {s} -> {result} ({reason})")