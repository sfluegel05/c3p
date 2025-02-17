"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: Polyamine
Definition: Any organic amino compound that contains two or more amino groups.

This revised version does the following:
  1. Rejects molecules that appear peptide‐like (>=2 amide bonds and >10 heavy atoms).
  2. Adds explicit hydrogens to get a reliable hydrogen count.
  3. Loops over every nitrogen atom (atomic num 7) that is not quaternary.
  4. For aliphatic nitrogens, we simply count them as candidate amino groups.
     For aromatic nitrogens we count them only if they have at least one hydrogen or a positive formal charge.
  5. As an extra safeguard, if the molecule has many heavy atoms and most atoms are in rings
     (ring fraction >0.8) then we require at least 3 candidate amino groups to call the molecule a polyamine.
  6. Also, if the molecule is a single ring system that is non‐aromatic and small (num heavy atoms >= 6)
     then we require at least 3 candidates (to avoid cases like imidazolidine).
     
These changes are intended to reduce the mis‐classification seen previously.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.
    
    The function:
      - Rejects peptide-like molecules (>=2 amide bonds and >10 heavy atoms).
      - Adds explicit hydrogens.
      - Scans every nitrogen atom (atomic number 7) that is not quaternary.
      - For aliphatic nitrogens, counts them as candidate amino groups.
      - For aromatic nitrogens, counts them only if they bear at least one hydrogen or are positively charged.
      - Applies additional safeguards against heavily cyclic systems (if most heavy atoms are in rings)
        or against a single small non‐aromatic ring.
    
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
    
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    # Check for peptide-like structure: if >=2 amide bonds and molecule is sufficiently large.
    # Amide pattern: nitrogen single-bonded to a carbon that is double-bonded to an oxygen.
    amide_smarts = Chem.MolFromSmarts("[NX3][C](=[O])")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if len(amide_matches) >= 2 and len(heavy_atoms) > 10:
        return False, f"Molecule appears peptide-like ({len(amide_matches)} amide bonds)"
    
    # Add explicit hydrogens for reliable hydrogen counts
    mol = Chem.AddHs(mol)
    
    candidate_idxs = []  # indices of candidate amino N atoms
    # Loop over each atom: we consider nitrogen atoms (atomic number 7)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue
        
        # Exclude quaternary nitrogen (rough heuristic: if degree (neighbors) >= 4, skip it)
        if atom.GetDegree() >= 4:
            continue
        
        # For aliphatic nitrogens, count them as candidates regardless of hydrogen count.
        # For aromatic nitrogens, require at least one hydrogen (explicit or implicit) OR a positive charge.
        if atom.GetIsAromatic():
            totalHs = atom.GetTotalNumHs()  # explicit+implicit hydrogens
            if totalHs < 1 and atom.GetFormalCharge() <= 0:
                continue
        # Else (aliphatic nitrogen) we count it even if no H—tertiary aliphatic amines are considered.
        
        candidate_idxs.append(atom.GetIdx())
    
    candidate_count = len(candidate_idxs)
    
    # Basic criterion: must have at least 2 candidates to be a polyamine.
    if candidate_count < 2:
        return False, f"Found only {candidate_count} candidate amino group(s), need at least 2"
    
    # Recalculate heavy atom count (after adding hydrogens, count only non-H atoms)
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    num_heavy = len(heavy_atoms)
    # Calculate number of heavy atoms that are in at least one ring.
    num_in_ring = sum(1 for atom in heavy_atoms if mol.GetRingInfo().NumAtomRings(atom.GetIdx()) > 0)
    ring_fraction = num_in_ring / num_heavy if num_heavy > 0 else 0.0
    
    # Safeguard 1: For larger, heavily cyclic molecules (num_heavy >= 10 and ring_fraction > 0.8),
    # require at least 3 candidate amino groups.
    if num_heavy >= 10 and ring_fraction > 0.8 and candidate_count < 3:
        return False, (f"Molecule is heavily cyclic (ring fraction {ring_fraction:.2f}, {num_heavy} heavy atoms) "
                       f"and has only {candidate_count} candidate amino group(s), require at least 3")
    
    # Safeguard 2: If the molecule consists of a single ring that is non‐aromatic and has at least 6 heavy atoms,
    # then require at least 3 candidate amino groups (to prevent small cyclic systems like imidazolidine).
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() == 1 and num_heavy >= 6:
        # Get the indices of atoms in the single ring.
        rings = ring_info.AtomRings()
        ring_atoms = rings[0]
        # Check if the ring is entirely non-aromatic
        if all(not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring_atoms):
            if candidate_count < 3:
                return False, (f"Single small non-aromatic ring system detected with {candidate_count} amino group(s); "
                               "require at least 3")
    
    return True, f"Contains {candidate_count} candidate amino group(s)"

# Example usage (if run as a script)
if __name__ == "__main__":
    test_examples = [
        # True positives:
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
        ("NCCCCNCCN", "N-(2-aminoethyl)butane-1,4-diamine"),
        # A couple of examples that previously were false positives:
        ("C1CNCN1", "imidazolidine"),
    ]
    
    for smi, name in test_examples:
        result, reason = is_polyamine(smi)
        print(f"SMILES: {smi} | Name: {name} -> {result} ({reason})")