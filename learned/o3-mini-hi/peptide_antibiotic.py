"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: Peptide antibiotics
A heuristic: a peptide antibiotic must be a valid molecule,
with a molecular weight of at least ~300 Da, contain several (≥3) amide bonds,
have a minimum nitrogen content (>7% of heavy atoms), and show at least one peptide‐like backbone fragment.
This is not perfect given the chemical diversity of peptide antibiotics.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is likely a peptide antibiotic based on its SMILES string.
    
    Heuristic criteria:
      - SMILES must parse to a valid molecule.
      - Molecular weight must be at least approximately 300 Da.
      - Must contain at least 3 amide bonds (defined here by a C(=O)N pattern).
      - Must have a minimum nitrogen content relative to heavy atoms (>7%).
      - Must contain at least one peptide-backbone fragment 
        (for example an alpha-carbon attached to an amine and a carbonyl: [C@H](N)C(=O) or a relaxed version if stereo is missing).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule likely belongs to the peptide antibiotic class, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a peptide antibiotic"
    
    # Look for amide bonds using a SMARTS pattern (matches C(=O)N).
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
    
    # Check for a peptide-backbone fragment.
    # First try with stereochemistry specified (common in peptides).
    peptide_pattern = Chem.MolFromSmarts("[C@H](N)C(=O)")
    peptide_matches = mol.GetSubstructMatches(peptide_pattern)
    # If not found, try a relaxed pattern without chiral information.
    if len(peptide_matches) < 1:
        peptide_pattern2 = Chem.MolFromSmarts("[C](N)C(=O)")
        peptide_matches = mol.GetSubstructMatches(peptide_pattern2)
    if len(peptide_matches) < 1:
        return False, "No typical peptide backbone fragment detected"
    
    # All criteria met, classify the molecule as a peptide antibiotic.
    return True, (f"Molecule weighs {mol_wt:.1f} Da with {n_amide} amide bonds, "
                  f"nitrogen ratio {ratio:.2f}, and peptide-like substructure present - "
                  "likely a peptide antibiotic")
    
# Example usage:
if __name__ == "__main__":
    # Test with surfactin C (should be classified as a peptide antibiotic)
    test_smiles = "[H][C@@]1(CCCCCCCCCC(C)C)CC(=O)N[C@@H](CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(=O)N[C@H](CC(C)C)C(=O)N[C@@H](C(C)C)C(=O)N[C@@H](CC(O)=O)C(=O)N[C@H](CC(C)C)C(=O)O1"
    result, reason = is_peptide_antibiotic(test_smiles)
    print(result, reason)