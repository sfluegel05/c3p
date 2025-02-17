"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: Peptide antibiotics
Heuristic improvements:
  - Molecule must have a molecular weight >= 300 Da.
  - Must have at least 3 amide bonds (using the substructure C(=O)N).
  - Nitrogen content among heavy atoms must be >7%.
  - For acyclic molecules:
      • If MW > 400 Da, require at least one contiguous dipeptide bond pattern (two consecutive amide bonds)
        detected by the relaxed SMARTS "[C](=O)[N]C(=O)[N]".
      • If MW is in the 300–400 Da range, require a longer contiguous sequence (three consecutive amide bonds)
        detected by the SMARTS "[C](=O)[N]C(=O)[N]C(=O)[N]".
  - If the molecule is cyclic, the contiguous pattern check is skipped.
This new heuristic aims to avoid misclassifying small (3‐residue) peptides while correctly recognizing 
peptide antibiotics, even when the backbone is embedded in a cyclic or complex topology.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is likely a peptide antibiotic based on its SMILES string.
    
    Heuristic criteria:
      - Valid molecule with molecular weight >= 300 Da.
      - At least 3 amide bonds (substructure "C(=O)N").
      - Nitrogen fraction (>7%) among heavy atoms.
      - For acyclic molecules:
          • Larger molecules (MW > 400 Da) must show at least one contiguous amide chain 
            (a dipeptide bond motif: "[C](=O)[N]C(=O)[N]").
          • Smaller molecules (300–400 Da) must show an even longer contiguous chain 
            (three consecutive amide bonds: "[C](=O)[N]C(=O)[N]C(=O)[N]") to reduce false positives.
      - For cyclic molecules, we skip the contiguous backbone check.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule likely belongs to the peptide antibiotic class, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate exact molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a peptide antibiotic"
    
    # Count amide bonds using a simple SMARTS pattern: C(=O)N
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    n_amide = len(amide_matches)
    if n_amide < 3:
        return False, f"Not enough amide bonds ({n_amide} found). Likely not a peptide chain"
    
    # Calculate nitrogen ratio among heavy atoms (non-hydrogen).
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    n_nitrogen = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 7)
    ratio = (n_nitrogen / len(heavy_atoms)) if heavy_atoms else 0
    if ratio < 0.07:
        return False, f"Low nitrogen ratio ({ratio:.2f}), not typical for peptides"
    
    # Determine whether the molecule is cyclic.
    is_cyclic = (mol.GetRingInfo().NumRings() > 0)
    
    # For acyclic molecules, require a contiguous chain of amide bonds.
    if not is_cyclic:
        if mol_wt > 400:
            # Look for two consecutive amide bonds (a dipeptide motif)
            contiguous_pattern = Chem.MolFromSmarts("[C](=O)[N]C(=O)[N]")
            matches = mol.GetSubstructMatches(contiguous_pattern)
            if len(matches) < 1:
                return False, "No contiguous dipeptide motif detected in a large acyclic molecule"
        else:
            # For molecules between 300-400 Da, be stricter to avoid small peptides.
            long_contiguous = Chem.MolFromSmarts("[C](=O)[N]C(=O)[N]C(=O)[N]")
            matches = mol.GetSubstructMatches(long_contiguous)
            if len(matches) < 1:
                return False, "No extended contiguous amide bond chain detected in a small acyclic molecule"
    # For cyclic molecules we assume that the multiple amide bonds and nitrogen content is evidence enough.
    
    # If all tests pass, return True with a summary.
    return True, (f"Molecule weighs {mol_wt:.1f} Da with {n_amide} amide bonds, nitrogen ratio {ratio:.2f}, "
                  f"{'cyclic' if is_cyclic else 'acyclic'} backbone detected - likely a peptide antibiotic")
    
# Example usage (for testing within a script):
if __name__ == "__main__":
    # Example: (+)-zwittermicin A: a known peptide antibiotic (expected True)
    test_smiles = "N[C@H](CO)[C@@H](O)C[C@H](O)[C@H](N)[C@H](O)[C@@H](O)C(=O)N[C@@H](CNC(N)=O)C(N)=O"
    result, reason = is_peptide_antibiotic(test_smiles)
    print(result, reason)