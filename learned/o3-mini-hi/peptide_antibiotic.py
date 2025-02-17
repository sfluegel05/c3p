"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: Peptide antibiotics

Improved heuristic criteria:
  • Molecule must be valid with molecular weight ≥ 300 Da.
  • Must have at least 3 amide bonds (found via SMARTS "C(=O)N").
  • The fraction of nitrogen (excluding hydrogens) must be at least 7%.
  • For cyclic molecules (i.e. at least one ring) the weight must be at least 400 Da.
  • Additionally, the molecule should contain at least one contiguous peptide-backbone motif,
    meaning that it includes a substructure like "[C@H](N)C(=O)" or its mirror.
    
This combination of criteria helps classify many peptide antibiotics while reducing false positives.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule likely belongs to the peptide antibiotic class based on its SMILES string.
    
    Heuristic criteria:
      - Valid molecule with molecular weight >= 300 Da.
      - Contains at least 3 amide bonds (SMARTS "C(=O)N").
      - Nitrogen fraction among heavy (non‐H) atoms is >= 0.07.
      - For cyclic molecules (has any ring), the molecular weight must be >= 400 Da.
      - In addition, at least one peptide-backbone motif (e.g. "[C@H](N)C(=O)" or "[C@@H](N)C(=O)") must be present.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is likely a peptide antibiotic, False otherwise.
        str: Explanation/reason for the classification decision.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate the exact molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a peptide antibiotic"
    
    # Count amide bonds (use a simple SMARTS pattern: C(=O)N ).
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    n_amide = len(amide_matches)
    if n_amide < 3:
        return False, f"Not enough amide bonds ({n_amide} found); typical peptide antibiotics have at least 3"
    
    # Calculate nitrogen ratio (excluding hydrogen atoms).
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    n_nitrogen = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 7)
    nitrogen_ratio = (n_nitrogen / len(heavy_atoms)) if heavy_atoms else 0
    if nitrogen_ratio < 0.07:
        return False, f"Low nitrogen ratio ({nitrogen_ratio:.2f}); not typical for peptides"
    
    # Check if the molecule is cyclic (has at least one ring).
    is_cyclic = (mol.GetRingInfo().NumRings() > 0)
    if is_cyclic and mol_wt < 400:
        return False, f"Cyclic molecule but weight too low ({mol_wt:.1f} Da) for a typical peptide antibiotic"
    
    # Now check for a peptide-backbone motif.
    # We look for a chiral carbon attached to an amino group and a carbonyl: "[C@H](N)C(=O)" or its mirror.
    pattern1 = Chem.MolFromSmarts("[C@H](N)C(=O)")
    pattern2 = Chem.MolFromSmarts("[C@@H](N)C(=O)")
    matches1 = mol.GetSubstructMatches(pattern1)
    matches2 = mol.GetSubstructMatches(pattern2)
    n_peptide_motifs = len(matches1) + len(matches2)
    if n_peptide_motifs < 1:
        return False, "No identifiable peptide backbone motif (e.g. '[C@H](N)C(=O)') found"
    
    backbone_type = "cyclic" if is_cyclic else "acyclic"
    explanation = (f"Molecule weighs {mol_wt:.1f} Da with {n_amide} amide bonds, "
                   f"nitrogen ratio {nitrogen_ratio:.2f}, {n_peptide_motifs} peptide motif(s), and a {backbone_type} backbone – "
                   "likely a peptide antibiotic.")
    return True, explanation

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Test with a known peptide antibiotic: (+)-zwittermicin A.
    test_smiles = "N[C@H](CO)[C@@H](O)C[C@H](O)[C@H](N)[C@H](O)[C@@H](O)C(=O)N[C@@H](CNC(N)=O)C(N)=O"
    is_pa, reason = is_peptide_antibiotic(test_smiles)
    print(is_pa, reason)