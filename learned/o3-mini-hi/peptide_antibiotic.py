"""
Classifies: CHEBI:25903 peptide antibiotic
"""
"""
Classifies: Peptide antibiotics
A heuristic: a peptide antibiotic must be a valid molecule
with a molecular weight of at least ~300 Da, contain several (≥3) amide bonds,
and have a significant fraction of nitrogen relative to the total heavy atoms.
Note: This approach will not be perfect given the chemical diversity of peptide antibiotics.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is likely a peptide antibiotic based on its SMILES string.
    
    The heuristic criteria used here:
      - SMILES must parse to a valid molecule.
      - The molecular weight must be at least 300 Da.
      - The molecule must contain at least 3 amide bonds (a proxy for peptide bonds).
      - The molecule must have a reasonable nitrogen content (>10% of heavy atoms),
        a typical feature for peptides.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule likely belongs to the peptide antibiotic class, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a peptide antibiotic"

    # Use a SMARTS pattern to detect amide bonds.
    # This pattern matches a carbonyl carbon attached to a nitrogen (C(=O)N).
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    n_amide = len(amide_matches)
    
    if n_amide < 3:
        return False, f"Not enough amide bonds ({n_amide} found). Likely not a peptide chain"

    # Calculate nitrogen content among heavy atoms (C, N, O, etc.).
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    n_nitrogen = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 7)
    
    if len(heavy_atoms) == 0 or (n_nitrogen / len(heavy_atoms)) < 0.10:
        ratio = (n_nitrogen / len(heavy_atoms)) if heavy_atoms else 0
        return False, f"Low nitrogen ratio ({ratio:.2f}), not typical for peptides"

    # If all criteria are met, we classify the molecule as a peptide antibiotic.
    return True, (f"Molecule weighs {mol_wt:.1f} Da with {n_amide} amide bonds "
                  f"and a nitrogen ratio of {n_nitrogen/len(heavy_atoms):.2f} - "
                  "likely a peptide antibiotic")