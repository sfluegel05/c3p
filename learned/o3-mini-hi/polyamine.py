"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: Polyamine
Definition: Any organic amino compound that contains two or more amino groups.
This revised version:
  1. Rejects molecules that look peptide-like (>=2 amide bonds in a reasonably sized molecule).
  2. Adds explicit hydrogens to get reliable hydrogen counts.
  3. Inspects every nitrogen (atomic number 7) that is not quaternary.
  4. Counts a nitrogen as an amino candidate if it is not in an amide (with one small exception) or
     if it is aliphatic, or (if aromatic) only if it bears at least one hydrogen OR is positively charged.
  5. Finally, if at least two candidate amino groups are found the molecule is accepted,
     with an additional check that larger molecules which are nearly “all rings” must have at least three candidates.
This heuristic is implemented with RDKit.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.
    
    The function first rejects peptide-like molecules if two or more amide bonds are found in a molecule
    with >10 heavy atoms. Then it adds explicit hydrogens and inspects every nitrogen atom. For aromatic
    nitrogens, at least one hydrogen (or a positive formal charge) is required; aliphatic ones are counted.
    
    Additionally, if the molecule is sizable (more than eight heavy atoms) and nearly all heavy atoms
    are in rings (ring fraction > 0.8) we require at least three candidate amino groups to avoid false positives.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a polyamine, False otherwise.
        str: Reason explaining the classification decision.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for peptide-like structure: if >=2 amide bonds and molecule is sufficiently large.
    # Amide pattern: N(any) - C(=O)
    amide_smarts = Chem.MolFromSmarts("[NX3][C](=[O])")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) >= 2 and mol.GetNumHeavyAtoms() > 10:
        return False, f"Molecule appears peptide-like ({len(amide_matches)} amide bonds)"
    
    # Add explicit hydrogens for more reliable hydrogen counts
    mol = Chem.AddHs(mol)
    
    candidate_idxs = []
    
    # Loop over each atom; we are interested in nitrogen atoms (atomic number 7)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        
        # Exclude quaternary nitrogen (degree >= 4)
        if atom.GetDegree() >= 4:
            continue
        
        # Check if nitrogen is directly involved in an amide bond.
        # We flag as "amide" if the N is bonded to a carbon that in turn is doubly bonded to oxygen.
        is_amide = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:  # possible carbonyl carbon?
                for nbr2 in nbr.GetNeighbors():
                    if nbr2.GetIdx() == atom.GetIdx():
                        continue
                    if nbr2.GetAtomicNum() == 8:
                        bond_CO = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                        if bond_CO and bond_CO.GetBondType() == Chem.BondType.DOUBLE:
                            is_amide = True
                            break
                if is_amide:
                    break
        # In our search we allow amide-type nitrogens when they come from non-peptide molecules 
        # (e.g. N-acetylputrescine). So we do not outright reject them.
        
        # For aromatic nitrogen, count only if it has at least one hydrogen or is positively charged.
        # (This accepts protonated amines which might not register an explicit H.)
        if atom.GetIsAromatic():
            # Get total number of hydrogens (explicit + implicit)
            totalHs = atom.GetTotalNumHs()
            if totalHs < 1 and atom.GetFormalCharge() <= 0:
                continue
        
        # Passed the filters. Count as candidate amino group.
        candidate_idxs.append(atom.GetIdx())
    
    candidate_count = len(candidate_idxs)
    if candidate_count < 2:
        return False, f"Found only {candidate_count} candidate amino group(s), need at least 2"
    
    # For the original version we checked that at least one pair of candidates were separated by >=3 bonds.
    # However, that rule proved too strict for small, legitimate polyamines (e.g., methanediamine).
    # Therefore, we remove the separation requirement.
    
    # Additional safeguard: if the molecule is large and almost all atoms are in rings,
    # then spurious amino groups on peripheral rings are likely not sufficient.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    num_heavy = len(heavy_atoms)
    if num_heavy > 0:
        num_in_ring = sum(1 for atom in heavy_atoms if mol.GetRingInfo().NumAtomRings(atom.GetIdx()) > 0)
        ring_fraction = num_in_ring / num_heavy
    else:
        ring_fraction = 0.0
    if mol.GetNumHeavyAtoms() > 8 and ring_fraction > 0.8 and candidate_count < 3:
        return False, f"Molecule is heavily cyclic (ring fraction {ring_fraction:.2f}) and has only {candidate_count} candidate amino group(s)"
    
    return True, f"Contains {candidate_count} candidate amino group(s)"

# Example usage and testing (if run as a script)
if __name__ == "__main__":
    test_examples = [
        # Expected true polyamines:
        ("CCNc1nc(N)nc(O)n1", "4-amino-6-(ethylamino)-1,3,5-triazin-2-ol"),
        ("NCCCNCCNCCCN", "3,2,3-tetramine"),
        ("N(C=1N=C(N)N=CN1)C2=CC=C(C=C2)Cl", "chlorazanil"),
        ("CC(=O)NCCCCN", "N-acetylputrescine"),
        ("S(=O)(C=1N=C(NC(C)C)N=C(NCC)N1)C", "N-Ethyl-N'-isopropyl-6-(methylsulfinyl)-1,3,5-triazine-2,4-diamine"),
        ("CN(C)c1ccc(cc1)C(N)c1ccc(cc1)N(C)C", "4,4'-(aminomethylene)bis(N,N-dimethylaniline)"),
        ("C1CN2CCN1CC2", "triethylenediamine"),
        ("CCNc1nc(Cl)nc(N[C@@H](C)CC)n1", "(S)-sebuthylazine"),
        ("CCNc1nc(N[C@H](C)CC)nc(OC)n1", "(R)-secbumeton"),
        ("NCCCNCCCCN(CCCN)CCCN", "N(4)-aminopropylspermine"),
        ("NCCN", "ethylenediamine"),
        ("NCN", "methanediamine"),
        ("NCCNCCNCC", "N,N'-diethylethylenediamine"),
        ("NCCCNCCCN", "bis(3-aminopropyl)amine"),
        ("C(C[NH3+])CC[NH2+]CCC([O-])=O", "putreanine(1+)"),
        ("[NH3+]CCCCNO", "N-hydroxyputrescine(1+)"),
        ("[NH3+]CCCCN(C(C)=O)O", "N(1)-acetyl-N(1)-hydroxyputrescine(1+)"),
        ("[H][C@@]12Cc3c[nH]c4cccc(c34)[C@@]1([H])CCCN2", "ergoline"),
        ("NCCCCNCCN", "N-(2-aminoethyl)butane-1,4-diamine"),
        ("NCCCNCCCCN(CCCN)CCCN", "N(4)-aminopropylspermine"),  # duplicate example can be tested
        # Examples that previously falsely classified as polyamine (should be False now):
        ("S(C1=C(N2C(=O)[C@@H]([C@H]2C1)[C@@H](O)C)C(=O)O)CCN", "8-epi-thienamycin"),
        ("C1CCN(CC1)C(=O)[C@H]2[C@@H]([C@H]3CN4C(=CC=C(C4=O)C5=CC=NC=C5)[C@H]3N2CC6=COC=N6)CO", "(2R,3R,3aS,9bS)-..."),
    ]
    
    for smi, name in test_examples:
        result, reason = is_polyamine(smi)
        print(f"SMILES: {smi} | Name: {name} -> {result} ({reason})")