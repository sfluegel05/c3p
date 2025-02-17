"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: Lipopolysaccharide
Definition:
    Lipopolysaccharides are natural compounds consisting of a trisaccharide repeating unit 
    (two heptose units and octulosonic acid) with oligosaccharide side chains and 3-hydroxytetradecanoic 
    acid units (major constituents of Gram-negative bacterial cell walls).
    
Heuristic criteria used:
  - At least 3 sugar-like rings (size 5–7 rings containing at least one oxygen and largely carbon atoms).
  - Presence of a carboxylic acid group (–C(=O)[OH]) present in sugar acids.
  - Presence of a long aliphatic chain (a chain of ≥12 CH2 groups) as a surrogate for fatty acid units.
  - Total molecular weight > 1000 Da.
  
Note: This is a heuristic classifier and may not capture all nuances of lipopolysaccharide structures.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as lipopolysaccharide, False otherwise.
        str: Reason explaining the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for sugar-like rings:
    # We count rings with size 5-7 that contain at least one oxygen and mostly carbon atoms.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_ring_count = 0
    for ring in ring_info:
        if len(ring) not in [5, 6, 7]:
            continue
        atom_symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
        # Count oxygen and carbon atoms.
        oxygen_count = atom_symbols.count("O")
        carbon_count = atom_symbols.count("C")
        # For a typical sugar ring, at least one oxygen must be present and most atoms should be carbons.
        if oxygen_count >= 1 and carbon_count >= (len(ring) - 1):
            sugar_ring_count += 1
    if sugar_ring_count < 3:
        return False, f"Found only {sugar_ring_count} sugar-like rings; at least 3 expected"

    # 2. Check for carboxylic acid groups.
    # SMARTS for carboxylic acid: a carbonyl carbon with a hydroxyl.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group (-C(=O)[OH]) found"

    # 3. Check for long aliphatic chain (as a proxy for fatty acid units).
    # We use a simple SMARTS that searches for at least 12 contiguous CH2 groups.
    fatty_chain_pattern = Chem.MolFromSmarts("[CH2]{12,}")
    if not mol.HasSubstructMatch(fatty_chain_pattern):
        return False, "No long aliphatic chain (≥12 CH2 groups) found to represent fatty acid unit"

    # 4. Check molecular weight (lipopolysaccharides are typically large, >1000 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for lipopolysaccharide"

    # If all heuristic criteria are met, classify as lipopolysaccharide.
    return True, ("Structure contains multiple sugar-like rings, carboxylic acid groups, "
                  "and long aliphatic chains consistent with a lipopolysaccharide")

# Example use (for testing purpose):
if __name__ == "__main__":
    # One of the provided SMILES examples for lipopolysaccharide (e.g., Mer-WF3010)
    test_smiles = "O=C(OC1C(O)C2(OC(C1OC3OC(C(O)C(C3O)O)COC(=O)/C=C/C=C/C=C/C)CO)OCC=4C2=C(O)C=C(O)C4)/C=C/C=C/CC(O)/C(=C/C=C/CCC(CC)C)/C"
    result, reason = is_lipopolysaccharide(test_smiles)
    print("Lipopolysaccharide classification:", result)
    print("Reason:", reason)