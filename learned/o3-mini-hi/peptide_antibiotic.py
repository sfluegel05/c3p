"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: Peptide antibiotics
A heuristic: A peptide antibiotic must be a valid molecule with a molecular weight 
of at least ~300 Da, contain several (≥3) amide bonds, have a minimal nitrogen content (>7%),
and possess a peptide‐like backbone. To reduce false positives, we now require detection 
of a contiguous dipeptide motif (“[C@H](N)C(=O)[C@H](N)”) in larger molecules (MW > 400 Da).
For smaller molecules (300–400 Da) a single amino‐acid fragment (“[C](N)C(=O)”) is allowed.
This improved heuristic better distinguishes true peptide antibiotics.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is likely a peptide antibiotic based on its SMILES string.
    
    Heuristic criteria:
      - The SMILES must produce a valid molecule.
      - Molecular weight must be at least approximately 300 Da.
      - The molecule must contain at least 3 amide bonds (C(=O)N substructure).
      - The nitrogen content (among heavy atoms) must be > 7%.
      - Larger molecules (MW > 400 Da) are required to have at least one 
        dipeptide motif ([C@H](N)C(=O)[C@H](N)) while smaller molecules must contain 
        a peptide-backbone fragment ([C](N)C(=O)).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule likely belongs to the peptide antibiotic class,
              False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a peptide antibiotic"
    
    # Count amide bonds using a SMARTS pattern "C(=O)N"
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    n_amide = len(amide_matches)
    if n_amide < 3:
        return False, f"Not enough amide bonds ({n_amide} found). Likely not a peptide chain"
    
    # Calculate nitrogen content among heavy (non-hydrogen) atoms.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    n_nitrogen = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 7)
    ratio = (n_nitrogen / len(heavy_atoms)) if heavy_atoms else 0
    if ratio < 0.07:
        return False, f"Low nitrogen ratio ({ratio:.2f}), not typical for peptides"
    
    # Check for peptide-backbone motifs.
    # For larger molecules (MW > 400 Da), require a contiguous dipeptide motif.
    dipeptide_pattern = Chem.MolFromSmarts("[C@H](N)C(=O)[C@H](N)")
    dipeptide_matches = mol.GetSubstructMatches(dipeptide_pattern)
    
    if mol_wt > 400:
        if len(dipeptide_matches) < 1:
            # If not found, try a more relaxed, non-stereospecific version.
            relaxed_dipeptide = Chem.MolFromSmarts("[C](N)C(=O)[C](N)")
            relaxed_matches = mol.GetSubstructMatches(relaxed_dipeptide)
            if len(relaxed_matches) < 1:
                return False, "No contiguous dipeptide motif detected in a large molecule"
    else:
        # For molecules near the weight cutoff, look for a single-unit peptide fragment.
        single_peptide = Chem.MolFromSmarts("[C](N)C(=O)")
        single_matches = mol.GetSubstructMatches(single_peptide)
        if len(single_matches) < 1:
            return False, "No peptide backbone fragment detected in a small molecule"
    
    # All criteria met; classify as a peptide antibiotic.
    return True, (f"Molecule weighs {mol_wt:.1f} Da with {n_amide} amide bonds, "
                  f"nitrogen ratio {ratio:.2f}, and peptide-like backbone detected - "
                  "likely a peptide antibiotic")
    
# Example usage if run as a script:
if __name__ == "__main__":
    # Test with surfactin C (expected True)
    test_smiles = "[H][C@@]1(CCCCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@H](CC(C)C)C(=O)O1"
    result, reason = is_peptide_antibiotic(test_smiles)
    print(result, reason)