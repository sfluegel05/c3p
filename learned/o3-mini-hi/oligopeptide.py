"""
Classifies: CHEBI:25676 oligopeptide
"""
"""
Classifies: Oligopeptide – a peptide containing a relatively small number of amino acids.
This improved heuristic uses several steps:
  - It searches for peptide (amide) bonds using the SMARTS "C(=O)N".
  - It estimates the number of amino acid residues as (number of peptide bonds + 1) and requires that
    count to be between 2 and 10.
  - It checks that sufficient alpha–carbons are present by counting occurrences of [C@H](N) or [C@@H](N).
    For a dipeptide exactly 2 such centers are required.
  - For peptides with 3 or more residues it enforces the presence of an internal backbone segment,
    using a relaxed SMARTS pattern "[CX3](N)C(=O)[CX3](N)" (i.e. without explicit chirality).
  - It ensures that all peptide–bond atoms are in a single contiguous fragment.
  - It checks that the overall molecular weight roughly matches the residue count (using 60–170 Da per residue).
Note: This heuristic does not capture every nuance of peptide chemistry.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oligopeptide(smiles: str):
    """
    Determines if a molecule is a small peptide (oligopeptide) from its SMILES string.
    The method applies several heuristics:
      - Searches for peptide (amide) bonds using a SMARTS pattern "C(=O)N".
      - Estimates the residue count as (# of peptide bonds + 1) and requires 2 to 10 residues.
      - Checks that enough alpha–carbon units (with chiral tag) are present.
          • For dipeptides exactly 2 alpha-carbons are needed.
      - For peptides with >=3 residues, requires an internal backbone segment
         using a pattern "[CX3](N)C(=O)[CX3](N)" (ignoring chiral tags).
      - Ensures that the atoms involved in the peptide bonds are in a single connected fragment.
      - Verifies that the molecular weight is roughly appropriate given the number of residues.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an oligopeptide, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find peptide (amide) bonds using a simple SMARTS pattern.
    peptide_bond_smarts = Chem.MolFromSmarts("C(=O)N")
    peptide_bond_matches = mol.GetSubstructMatches(peptide_bond_smarts)
    n_peptide_bonds = len(peptide_bond_matches)
    
    if n_peptide_bonds == 0:
        return False, "No peptide (amide) bonds found"

    # Estimate the number of residues.
    n_residues = n_peptide_bonds + 1
    if n_residues < 2:
        return False, "Too few peptide bonds to form a peptide (need at least 2 residues)"
    if n_residues > 10:
        return False, f"Found {n_residues} amino acid residues which is too many for an oligopeptide"
    
    # Count occurrences of alpha–carbon atoms with chiral tags.
    alpha_smarts1 = Chem.MolFromSmarts("[C@H](N)")
    alpha_smarts2 = Chem.MolFromSmarts("[C@@H](N)")
    alpha_matches1 = mol.GetSubstructMatches(alpha_smarts1)
    alpha_matches2 = mol.GetSubstructMatches(alpha_smarts2)
    alpha_atoms = set()
    for match in alpha_matches1:
        alpha_atoms.add(match[0])
    for match in alpha_matches2:
        alpha_atoms.add(match[0])
    n_alpha = len(alpha_atoms)
    
    # For a dipeptide, require exactly 2 alpha–carbons.
    if n_residues == 2 and n_alpha != 2:
        return False, f"For a dipeptide, exactly 2 alpha–carbon centers are expected; found {n_alpha}"
    # For longer peptides, require at least 2 alpha–carbons.
    if n_residues >= 3 and n_alpha < 2:
        return False, f"Expected at least 2 alpha–carbon centers; found {n_alpha}"

    # For peptides with 3 or more residues, require an internal backbone segment.
    # Use a more relaxed SMARTS (ignoring chirality) that should match consecutive amino acids.
    if n_residues >= 3:
        internal_backbone = Chem.MolFromSmarts("[CX3](N)C(=O)[CX3](N)")
        if not mol.HasSubstructMatch(internal_backbone):
            return False, "Peptide backbone internal segment pattern not found"
    
    # Ensure that all atoms involved in the peptide bonds (the carbonyl carbon and the adjacent nitrogen)
    # belong to a single contiguous fragment.
    pb_atoms = set()
    for match in peptide_bond_matches:
        # The SMARTS "C(=O)N" returns three atoms:
        # match[0] = carbonyl carbon, match[1] = oxygen, match[2] = neighboring nitrogen.
        pb_atoms.add(match[0])
        pb_atoms.add(match[2])
    
    # Get the connected fragments (each fragment is a tuple of atom indices).
    fragments = Chem.GetMolFrags(mol, asMols=False)
    n_frag_with_pb = sum(1 for frag in fragments if any(atom in frag for atom in pb_atoms))
    if n_frag_with_pb > 1:
        return False, "Peptide bonds are not contiguous in a single backbone"
    
    # Calculate additional properties.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Set rough molecular weight bounds per residue (60–170 Da per residue).
    min_expected = n_residues * 60
    max_expected = n_residues * 170
    if mol_wt < min_expected:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for {n_residues} residue peptide (expected at least {min_expected} Da)"
    if mol_wt > max_expected:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too high for {n_residues} residue peptide (expected at most {max_expected} Da)"
    
    reason = (f"Detected {n_peptide_bonds} peptide bond(s) (≈{n_residues} residue(s)), {n_rotatable} rotatable bonds, "
              f"{n_alpha} alpha–carbon center(s), and MW of {mol_wt:.1f} Da. "
              "Peptide-bond atoms are in a single contiguous fragment and an appropriate backbone segment was identified. "
              "This is consistent with an oligopeptide.")
    
    return True, reason