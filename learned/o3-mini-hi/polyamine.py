"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: Polyamine
Definition: Any organic amino compound that contains two or more amino groups.
For our purposes we (a) try to rule out peptide‐like compounds (with multiple amide bonds),
(b) note that even if one nitrogen is “acetylated” sometimes that amine may still be counted in non‐peptide like molecules,
(c) count as an amino group any nitrogen that is not “locked up” in a fully aromatic heterocycle and that is not quaternary.
We also require that at least one pair of candidate nitrogen atoms is connected by a bond‐path of at least 3
(beyond very small molecules) to help rule out very compact cyclic systems.
This is a heuristic approach using RDKit.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.
    (An amino group is here defined in a heuristic way.)
    
    The function first flags molecules that look peptide-like (2 or more amide bonds in a >10 heavy atom molecule)
    and returns False in that case.
    Then it adds explicit hydrogens and inspects every nitrogen.
    
    For each nitrogen:
      - If it is directly attached (via a double bond) to a carbonyl carbon we call it an “amide” nitrogen.
        If the molecule looks peptide‐like then such N’s are not counted.
      - If the nitrogen is non‐ring, it is counted.
      - If the nitrogen is in a ring:
            * if it is fully aromatic, then we only count it if it bears at least one hydrogen (eg an –NH group)
              (pyridinic N without H are not counted).
            * if it is non‐aromatic, we count it, but later we check the connectivity of candidate amino groups.
    Finally, if at least two candidate amino groups were found,
    we require that at least one pair is “sufficiently separated” (path length >= 3, using all bonds) – if the molecule
    is acyclic. For molecules where every heavy atom is in a ring, we require that the heavy atom count be greater than 6.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a polyamine, False otherwise.
        str: A reason explaining the result.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for peptide-like: at least 2 amide bonds in a reasonably sized molecule.
    amide_smarts = Chem.MolFromSmarts("[NX3][C](=[O])")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) >= 2 and mol.GetNumHeavyAtoms() > 10:
        return False, f"Molecule appears peptide-like ({len(amide_matches)} amide bonds)"
    
    # Add explicit hydrogens so that we can reliably count attached H's
    mol = Chem.AddHs(mol)
    
    candidate_idxs = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # not nitrogen

        # Exclude quaternary N (degree >= 4) because these are ammonium centers
        if atom.GetDegree() >= 4:
            continue

        # Determine if this nitrogen is directly involved in an amide bond.
        # (If it is, and the molecule looked peptide-like then we already returned False.
        #  Otherwise, in a non-peptide molecule we count it anyway)
        is_amide = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                for bond in mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()):
                    # Look for double bond to oxygen on the carbon neighbor
                    # (skip if the oxygen is the original atom)
                    for c_nbr in nbr.GetNeighbors():
                        if c_nbr.GetAtomicNum() == 8 and c_nbr.GetIdx() != atom.GetIdx():
                            if mol.GetBondBetweenAtoms(nbr.GetIdx(), c_nbr.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                                is_amide = True
                                break
                    if is_amide:
                        break
            if is_amide:
                break
        # In our non-peptide molecules we do count an amide nitrogen.
        # (This allows e.g. N-acetylputrescine to still be a polyamine.)
        
        # Now decide by ring/aromatic environment:
        if atom.GetIsAromatic():
            # For aromatic N, count only if it has one or more attached hydrogens.
            # (This will skip pyridine N but include e.g. an indole NH.)
            if atom.GetTotalNumHs() < 1:
                continue  # skip aromatic nitrogen lacking any H
            # Otherwise, count it.
        # For non-aromatic ones, count them (whether they bear explicit H or are tertiary amines)
        candidate_idxs.append(atom.GetIdx())
    
    candidate_count = len(candidate_idxs)
    if candidate_count < 2:
        return False, f"Found only {candidate_count} candidate amino group(s), need at least 2"
    
    # Next, check if the recorded candidate amine groups are separated by at least two intervening bonds.
    # We compute the pairwise shortest path lengths (using all bonds) between candidate Ns.
    # For acyclic molecules (or molecules with some non‐ring bonds connecting the amine sites)
    # we require at least one pair to have a path length of >= 3.
    # For molecules where every heavy atom is in a ring, we require that the molecule is not extremely small.
    ring_info = mol.GetRingInfo()
    all_in_ring = all([ring_info.NumAtomRings(a.GetIdx()) > 0 for a in mol.GetAtoms() if a.GetAtomicNum() > 1])
    separation_ok = False
    # Use the full (cyclic) bond graph for distances.
    dmat = Chem.GetDistanceMatrix(mol)
    for i in range(len(candidate_idxs)):
        for j in range(i+1, len(candidate_idxs)):
            d = dmat[candidate_idxs[i]][candidate_idxs[j]]
            if d >= 3:
                separation_ok = True
                break
        if separation_ok:
            break
    if not separation_ok:
        if all_in_ring and mol.GetNumHeavyAtoms() > 6:
            # In a larger fused cyclic system we can sometimes have close candidate N’s
            separation_ok = True
        else:
            return False, "Candidate amino groups are too close together (likely in a compact cyclic system)"
    
    return True, f"Contains {candidate_count} candidate amino groups and at least one pair separated by >=3 bonds."

# Example usage and testing
if __name__ == "__main__":
    test_examples = [
        # True positives (expected polyamines)
        ("CCNc1nc(N)nc(O)n1", "4-amino-6-(ethylamino)-1,3,5-triazin-2-ol"),
        ("NCCCNCCNCCCN", "3,2,3-tetramine"),
        ("N(C=1N=C(N)N=CN1)C2=CC=C(C=C2)Cl", "chlorazanil"),
        ("S(=O)(C=1N=C(NC(C)C)N=C(NCC)N1)C", "N-Ethyl-N'-isopropyl-6-(methylsulfinyl)-1,3,5-triazine-2,4-diamine"),
        ("CN(C)c1ccc(cc1)C(N)c1ccc(cc1)N(C)C", "4,4'-(aminomethylene)bis(N,N-dimethylaniline)"),
        ("C1CN2CCN1CC2", "triethylenediamine"),
        ("CCNc1nc(Cl)nc(N[C@@H](C)CC)n1", "(S)-sebuthylazine"),
        ("CCNc1nc(N[C@H](C)CC)nc(OC)n1", "(R)-secbumeton"),
        ("NCCCNCCCCN(CCCN)CCCN", "N(4)-aminopropylspermine"),
        ("C(C[NH3+])CC[NH2+]CCC([O-])=O", "putreanine(1+)"),
        ("NCCN", "ethylenediamine"),
        # (Additional examples from the provided list could be added.)
        # False positives (should not be classified as polyamine):
        ("C1CNCN1", "imidazolidine"),  # small saturated cycle – too compact
    ]
    
    for smi, name in test_examples:
        result, reason = is_polyamine(smi)
        print(f"SMILES: {smi} | Name: {name} -> {result} ({reason})")