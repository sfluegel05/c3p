"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: Polyamine
Definition: Any organic amino compound that contains two or more amino groups.
A polyamine is here detected by:
  1. Rejecting molecules that look peptide-like (>=2 amide bonds in a >10 heavy atom molecule).
  2. Adding explicit hydrogens and then inspecting every nitrogen.
  3. For each nitrogen, excluding quaternary centers and (for aromatic nitrogens) those lacking an –NH.
  4. Also flagging nitrogen involved in an amide bond (if attached to a carbon that is double-bonded to oxygen).
  5. Finally, checking that at least two candidate amino groups exist and that at least one pair is sufficiently separated in the bond network.
This heuristic uses RDKit.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.
    
    The function first flags molecules that look peptide-like (>=2 amide bonds in a reasonably sized molecule)
    and returns False if so. Then it adds explicit hydrogens and inspects every nitrogen.
    
    For each nitrogen:
      - If it is directly attached (via a double bond) to a carbonyl carbon, we flag it as being in an amide.
      - Exclude quaternary nitrogens (degree >=4).
      - For aromatic nitrogen, count only if it has one or more hydrogens attached.
      - For non-aromatic nitrogens, simply count them.
    
    Finally, if at least two candidate amino groups are found, the code checks if at least one pair is separated
    by a bond path of at least 3 (using the full bond graph). For molecules in which every heavy atom is in a ring,
    a less strict requirement is applied provided the molecule is not extremely small.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a polyamine, False otherwise.
        str: A reason explaining the result.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for peptide-like structure: if >=2 amide bonds and molecule is sufficiently large.
    amide_smarts = Chem.MolFromSmarts("[NX3][C](=[O])")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) >= 2 and mol.GetNumHeavyAtoms() > 10:
        return False, f"Molecule appears peptide-like ({len(amide_matches)} amide bonds)"
    
    # Add explicit hydrogens for more reliable hydrogen counts
    mol = Chem.AddHs(mol)
    
    candidate_idxs = []
    
    # Traverse each atom checking for nitrogen candidates.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # Skip non-nitrogen atoms
        
        # Exclude quaternary nitrogen (degree >= 4) as these are ammonium centers.
        if atom.GetDegree() >= 4:
            continue
        
        # Check if the nitrogen is directly involved in an amide bond.
        is_amide = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # Neighbor is a carbon
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                # For the carbon neighbor, check if any of its other neighbors is an oxygen
                # involved in a double bond.
                for carbon_nbr in nbr.GetNeighbors():
                    if carbon_nbr.GetIdx() == atom.GetIdx():
                        continue  # Skip the original nitrogen atom
                    if carbon_nbr.GetAtomicNum() == 8:
                        bond_CO = mol.GetBondBetweenAtoms(nbr.GetIdx(), carbon_nbr.GetIdx())
                        if bond_CO and bond_CO.GetBondType() == Chem.BondType.DOUBLE:
                            is_amide = True
                            break
                if is_amide:
                    break
        
        # In our non-peptide molecules we count even amide nitrogens (e.g., N-acetylputrescine).
        # Decide based on aromaticity:
        if atom.GetIsAromatic():
            # For aromatic N, count only if it has at least one attached hydrogen (e.g. -NH)
            if atom.GetTotalNumHs() < 1:
                continue
        # If the nitrogen passes these checks, record its index.
        candidate_idxs.append(atom.GetIdx())
    
    candidate_count = len(candidate_idxs)
    if candidate_count < 2:
        return False, f"Found only {candidate_count} candidate amino group(s), need at least 2"
    
    # Next, ensure that at least one pair of candidate nitrogens is sufficiently separated.
    # This is computed by using the full (cyclic) bond distance matrix.
    ring_info = mol.GetRingInfo()
    # Identify if every heavy atom is in a ring.
    all_in_ring = all([ring_info.NumAtomRings(a.GetIdx()) > 0 for a in mol.GetAtoms() if a.GetAtomicNum() > 1])
    
    separation_ok = False
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
            separation_ok = True  # Allow in larger fused cyclic systems.
        else:
            return False, "Candidate amino groups are too close together (likely in a compact cyclic system)"
    
    return True, f"Contains {candidate_count} candidate amino groups and at least one pair separated by >=3 bonds."

# Example usage and testing (if run as a script)
if __name__ == "__main__":
    test_examples = [
        # Examples expected to be classified as polyamines:
        ("CCNc1nc(N)nc(O)n1", "4-amino-6-(ethylamino)-1,3,5-triazin-2-ol"),
        ("NCCCNCCNCCCN", "3,2,3-tetramine"),
        ("N(C=1N=C(N)N=CN1)C2=CC=C(C=C2)Cl", "chlorazanil"),
        ("CC(=O)NCCCCN", "N-acetylputrescine"),
        ("CN(C)c1ccc(cc1)C(N)c1ccc(cc1)N(C)C", "4,4'-(aminomethylene)bis(N,N-dimethylaniline)"),
        ("C1CN2CCN1CC2", "triethylenediamine"),
        ("CCNc1nc(Cl)nc(N[C@@H](C)CC)n1", "(S)-sebuthylazine"),
        ("CCNc1nc(N[C@H](C)CC)nc(OC)n1", "(R)-secbumeton"),
        ("NCCCNCCCCN(CCCN)CCCN", "N(4)-aminopropylspermine"),
        ("NCCN", "ethylenediamine"),
        ("NCN", "methanediamine"),
        ("NCCNCCNCC", "N,N'-diethylethylenediamine"),
        ("NCCCNCCCN", "bis(3-aminopropyl)amine"),
        # Example not expected to be a polyamine (compact cyclic system)
        ("C1CNCN1", "imidazolidine"),
    ]
    
    for smi, name in test_examples:
        result, reason = is_polyamine(smi)
        print(f"SMILES: {smi} | Name: {name} -> {result} ({reason})")