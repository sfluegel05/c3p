"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: Lipopolysaccharide
Definition:
    Lipopolysaccharides are natural compounds consisting of a trisaccharide repeating unit 
    (two heptose units and octulosonic acid) with oligosaccharide side chains and 3-hydroxytetradecanoic 
    acid units (major constituents of Gram-negative bacterial cell walls).

Heuristic criteria used in this classifier (improved):
  - At least one sugar-like ring is detected (rings of size 5–7 atoms containing ≥1 oxygen and mostly carbons).
    (Note: ideally a full repeating trisaccharide would have 3 sugars, but many examples presented fewer clear rings.)
  - Presence of a carboxyl/carboxylic acid motif is required. To capture both free acids and esterified forms,
    we look for a carbonyl carbon attached to an oxygen ([CX3](=O)[OX2]).
  - Presence of a long aliphatic chain. Here we search for at least 8 contiguous CH2 groups (relaxed from 12)
    to account for minor interruptions (e.g. hydroxylation) in 3-hydroxytetradecanoic acid units.
  - Molecular weight is expected to be >1000 Da.
    
Note:
  This is a heuristic approach; some lipopolysaccharides may deviate from these criteria.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a lipopolysaccharide, False otherwise.
        str: Reason explaining the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for sugar-like rings.
    # We count rings of size 5-7 that have at least one oxygen and mostly carbon atoms.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_ring_count = 0
    for ring in ring_info:
        if len(ring) not in [5, 6, 7]:
            continue
        # get the element symbols for atoms in the ring
        symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
        oxygen_count = symbols.count("O")
        carbon_count = symbols.count("C")
        # For a typical sugar ring, expect at least one oxygen and nearly all atoms to be carbons.
        if oxygen_count >= 1 and carbon_count >= (len(ring) - 1):
            sugar_ring_count += 1
    # Relaxed threshold: if at least one ring is found, we pass this criterion.
    if sugar_ring_count < 1:
        return False, f"Found only {sugar_ring_count} sugar-like rings; at least 1 expected"
    
    # 2. Check for carboxyl or ester-like moieties.
    # SMARTS pattern that matches a carbonyl carbon attached to an oxygen (matches both acid and ester forms).
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxyl/carboxylic acid motif ([CX3](=O)[OX2]) found"
    
    # 3. Look for a long aliphatic chain.
    # We search for at least 8 contiguous -CH2- groups (relaxed from 12 to allow for common modifications).
    fatty_chain_pattern = Chem.MolFromSmarts("[CH2]{8,}")
    if not mol.HasSubstructMatch(fatty_chain_pattern):
        return False, "No long aliphatic chain (≥8 contiguous CH2 groups) found to represent fatty acid units"
    
    # 4. Check the molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a lipopolysaccharide"

    return True, ("Structure contains sugar-like ring(s), carboxyl/ester motifs, "
                  "and long aliphatic chains with a sufficiently high molecular weight, "
                  "consistent with a lipopolysaccharide. (Note: a full trisaccharide repeating unit is ideal)")
    
# Example use (for testing purposes):
if __name__ == "__main__":
    # Example SMILES from one of the provided lipopolysaccharide examples (Mer-WF3010)
    test_smiles = "O=C(OC1C(O)C2(OC(C1OC3OC(C(O)C(C3O)O)COC(=O)/C=C/C=C/C=C/C)CO)OCC=4C2=C(O)C=C(O)C4)/C=C/C=C/CC(O)/C(=C/C=C/CCC(CC)C)/C"
    result, reason = is_lipopolysaccharide(test_smiles)
    print("Lipopolysaccharide classification:", result)
    print("Reason:", reason)