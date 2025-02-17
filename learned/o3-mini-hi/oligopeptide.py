"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: Oligopeptide – a peptide containing a relatively small number of amino acids.
This improved heuristic function identifies peptide (amide) bonds but now also
(+1) insists that these bonds occur as one contiguous backbone,
(+2) checks for at least one alpha–carbon unit (using chiral tags),
(+3) restricts the overall “residue count” (peptide bonds + 1) to no more than about 10,
and (+4) roughly checks that the molecular weight is in an expected range for that many residues.
Note: This heuristic will not capture every nuance of peptide chemistry.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is an oligopeptide from its SMILES string.
    The method uses several heuristics:
      - It searches for peptide (amide) bonds using the SMARTS "C(=O)N".
      - It estimates the residue count as (number of peptide bonds + 1)
        and requires that count to be between 2 and 10.
      - It verifies that at least one typical alpha–carbon unit (with chiral tag)
        is present (looking for [C@H](N) or [C@@H](N)).
      - It checks that all peptide-bond atoms lie in one connected fragment
        (to ensure the peptide backbone is contiguous).
      - It also checks that the molecular weight is in a plausible range
        (roughly between 70 and 150 Da per residue).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an oligopeptide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Search for amide (peptide) bonds.
    # Note: The pattern "C(=O)N" is intentionally simple; it will match many amide bonds.
    peptide_bond_smarts = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_smarts)
    n_peptide_bonds = len(peptide_bond_matches)
    
    if n_peptide_bonds == 0:
        return False, "No peptide (amide) bonds found"
    
    # Estimate number of amino acid residues (peptide bonds + 1).
    n_residues = n_peptide_bonds + 1
    if n_residues < 2:
        return False, "Too few peptide bonds to form a peptide (need at least 2 residues)"
    if n_residues > 10:
        return False, f"Found {n_residues} amino acid residues which is too many for an oligopeptide"
    
    # Check for the presence of at least one typical alpha–carbon unit.
    aa_pattern1 = Chem.MolFromSmarts("[C@H](N)")   # chiral with NH attached
    aa_pattern2 = Chem.MolFromSmarts("[C@@H](N)")
    has_alpha = mol.HasSubstructMatch(aa_pattern1) or mol.HasSubstructMatch(aa_pattern2)
    if not has_alpha:
        return False, "No typical alpha–amino acid unit (with chiral center) found"
    
    # Verify that all atoms involved in peptide bonds occur in the same connected fragment.
    # (We build a set of atom indices from all peptide bond matches.)
    pb_atoms = set()
    for match in peptide_bond_matches:
        # each match returns a tuple of (carbonyl carbon, oxygen, nitrogen) for "C(=O)N"
        # We are interested in the carbon and the nitrogen as they connect residues.
        pb_atoms.add(match[0])
        pb_atoms.add(match[2])
    
    # Get fragment information (each fragment is a tuple of atom indices).
    frags = Chem.GetMolFrags(mol, asMols=False)
    # See in how many fragments the peptide-bond atoms fall.
    frag_hits = [any(atom in frag for atom in pb_atoms) for frag in frags]
    # Count fragments that include at least one peptide-bond atom.
    n_frag_with_pb = sum(1 for frag in frags if any(atom in frag for atom in pb_atoms))
    if n_frag_with_pb > 1:
        return False, "Peptide bonds are not contiguous in a single backbone"
    
    # Calculate additional properties.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Make a rough check: typical amino acid residues are around 70–150 Da.
    min_expected = n_residues * 70
    max_expected = n_residues * 150
    if mol_wt < min_expected:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for {n_residues} residues (expected at least {min_expected} Da)"
    if mol_wt > max_expected:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too high for {n_residues} residues (expected at most {max_expected} Da)"
    
    reason = (f"Detected {n_peptide_bonds} peptide bonds (≈{n_residues} residues), "
              f"{n_rotatable} rotatable bonds, and MW of {mol_wt:.1f} Da. "
              "All peptide-bond atoms occur in one contiguous fragment and the weight is in the expected range. "
              "This is consistent with an oligopeptide.")
    
    return True, reason