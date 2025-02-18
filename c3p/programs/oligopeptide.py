"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: Oligopeptide – a peptide containing a relatively small number of amino acids.
This heuristic identifies peptide (amide) bonds, checks that these occur in one contiguous backbone,
verifies the presence of at least one typical alpha–carbon unit (with chiral tag),
ensures that, if more than two residues are present, an internal backbone segment is found,
and also roughly checks that the molecular weight is within a reasonable range.
Note: This heuristic does not capture every nuance of peptide chemistry.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is a small peptide (oligopeptide) from its SMILES string.
    The method uses several heuristics:
      - It searches for peptide (amide) bonds using the SMARTS "C(=O)N".
      - It estimates the residue count as (number of peptide bonds + 1) and requires that
        count to be between 2 and 10.
      - It verifies that at least one typical alpha–carbon (with chiral tag) is present,
        using SMARTS patterns [C@H](N) and [C@@H](N).
      - For peptides with ≥3 residues it enforces that an internal backbone segment is present,
        via a search for a pattern of two consecutively linked residues (e.g. "[C@H](N)C(=O)[C@H](N)")
      - It checks that all peptide–bond atoms fall in a single connected fragment.
      - It also checks that the overall molecular weight roughly matches the residue
        count using thresholds (approximately 60–170 Da per residue).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an oligopeptide, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify amide (peptide) bonds; note that "C(=O)N" is a simple pattern that will match many amide bonds.
    peptide_bond_smarts = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_smarts)
    n_peptide_bonds = len(peptide_bond_matches)
    
    if n_peptide_bonds == 0:
        return False, "No peptide (amide) bonds found"
    
    # Estimate the number of amino acid residues.
    n_residues = n_peptide_bonds + 1
    if n_residues < 2:
        return False, "Too few peptide bonds to form a peptide (need at least 2 residues)"
    if n_residues > 10:
        return False, f"Found {n_residues} amino acid residues which is too many for an oligopeptide"
    
    # Check for the presence of at least one typical alpha–carbon unit (with a chiral tag)
    aa_pattern1 = Chem.MolFromSmarts("[C@H](N)")
    aa_pattern2 = Chem.MolFromSmarts("[C@@H](N)")
    if not (mol.HasSubstructMatch(aa_pattern1) or mol.HasSubstructMatch(aa_pattern2)):
        return False, "No typical alpha–amino acid unit (with chiral center) found"
    
    # For peptides with three or more residues, require that an internal backbone segment is present.
    # This checks for two consecutive amino acid units connected by a peptide bond.
    if n_residues >= 3:
        internal_patterns = [
            Chem.MolFromSmarts("[C@H](N)C(=O)[C@H](N)"),
            Chem.MolFromSmarts("[C@H](N)C(=O)[C@@H](N)"),
            Chem.MolFromSmarts("[C@@H](N)C(=O)[C@H](N)"),
            Chem.MolFromSmarts("[C@@H](N)C(=O)[C@@H](N)")
        ]
        internal_match = False
        for patt in internal_patterns:
            if mol.HasSubstructMatch(patt):
                internal_match = True
                break
        if not internal_match:
            return False, "Peptide backbone internal segment pattern not found"
    
    # Ensure that all atoms involved in the peptide bonds (the carbonyl carbon and the adjacent nitrogen)
    # are located in a single connected fragment.
    pb_atoms = set()
    for match in peptide_bond_matches:
        # Each match for "C(=O)N" returns a tuple. We take the carbonyl carbon (match[0])
        # and the nitrogen (match[2]); these are the connecting atoms between residues.
        pb_atoms.add(match[0])
        pb_atoms.add(match[2])
    
    # Get the molecule fragmentations.
    frags = Chem.GetMolFrags(mol, asMols=False)
    n_frag_with_pb = sum(1 for frag in frags if any(atom in frag for atom in pb_atoms))
    if n_frag_with_pb > 1:
        return False, "Peptide bonds are not contiguous in a single backbone"
    
    # Calculate additional properties.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Set rough molecular weight bounds per residue.
    # (The lower bound is set to 60 Da/residue and the upper bound is 170 Da/residue.)
    min_expected = n_residues * 60
    max_expected = n_residues * 170
    if mol_wt < min_expected:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for {n_residues} residues (expected at least {min_expected} Da)"
    if mol_wt > max_expected:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too high for {n_residues} residues (expected at most {max_expected} Da)"
    
    reason = (f"Detected {n_peptide_bonds} peptide bonds (≈{n_residues} residues), "
              f"{n_rotatable} rotatable bonds, and MW of {mol_wt:.1f} Da. "
              "All peptide-bond atoms occur in one contiguous fragment and an internal backbone segment was identified. "
              "The molecular weight is in the expected range. This is consistent with an oligopeptide.")
    
    return True, reason