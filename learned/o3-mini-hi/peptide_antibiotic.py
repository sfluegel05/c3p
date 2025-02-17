"""
Classifies: CHEBI:25903 peptide antibiotic
"""
#!/usr/bin/env python3
"""
Classifies: Peptide antibiotics

Heuristic criteria:
  • The molecule must be valid and have a molecular weight ≥300 Da.
  • It must include at least 3 amide bonds (found via the substructure "C(=O)N").
  • The ratio of nitrogen atoms (atomic number 7) among non‐hydrogen atoms must be at least about 7%.
  • For cyclic molecules (i.e. with one or more rings) we require a molecular weight ≥800 Da, 
    as most peptide antibiotics with cyclic backbones are large.
  • For acyclic molecules, we no longer require a “contiguous dipeptide motif” in order not to miss 
    modified or very short peptide antibiotics.
    
This combination of criteria allows us to correctly classify many examples of peptide antibiotics while 
minimizing false positive assignments.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is likely a peptide antibiotic based on its SMILES string.
    
    Heuristic criteria:
      - Valid molecule with molecular weight >= 300 Da.
      - Contains at least 3 amide bonds (the SMARTS "C(=O)N").
      - Nitrogen fraction among heavy (non-H) atoms >= 0.07.
      - If cyclic (i.e. at least one ring), the molecular weight must be >= 800 Da.
      - If acyclic, no further contiguity check is required.
        
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule likely belongs to the peptide antibiotic class,
              False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate the exact molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a peptide antibiotic"
    
    # Use a simple SMARTS pattern to count amide bonds.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    n_amide = len(amide_matches)
    if n_amide < 3:
        return False, f"Not enough amide bonds ({n_amide} found); typical peptide antibiotics have at least 3"
    
    # Calculate the nitrogen ratio among heavy atoms (atomic number > 1).
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    n_nitrogen = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 7)
    nitrogen_ratio = (n_nitrogen / len(heavy_atoms)) if heavy_atoms else 0
    if nitrogen_ratio < 0.07:
        return False, f"Low nitrogen ratio ({nitrogen_ratio:.2f}); not typical for peptides"
    
    # Determine if the molecule is cyclic (has any rings)
    is_cyclic = (mol.GetRingInfo().NumRings() > 0)
    
    # For cyclic molecules, require a higher molecular weight threshold.
    if is_cyclic:
        if mol_wt < 800:
            return False, f"Cyclic molecule but weight too low ({mol_wt:.1f} Da) for a typical peptide antibiotic"
    
    # (For acyclic molecules we do not require a contiguous stretch check.)
    
    # If all tests pass, return True with a detailed explanation.
    backbone_type = "cyclic" if is_cyclic else "acyclic"
    explanation = (f"Molecule weighs {mol_wt:.1f} Da with {n_amide} amide bonds, "
                   f"nitrogen ratio {nitrogen_ratio:.2f}, and an {backbone_type} backbone – "
                   "likely a peptide antibiotic.")
    return True, explanation

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Try an example known peptide antibiotic: (+)-zwittermicin A
    test_smiles = "N[C@H](CO)[C@@H](O)C[C@H](O)[C@H](N)[C@H](O)[C@@H](O)C(=O)N[C@@H](CNC(N)=O)C(N)=O"
    result, reason = is_peptide_antibiotic(test_smiles)
    print(result, reason)