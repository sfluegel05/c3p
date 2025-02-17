"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: Lipopolysaccharide

Definition:
    Lipopolysaccharides are natural compounds consisting of a trisaccharide repeating unit 
    (two heptose units and octulosonic acid) with oligosaccharide side chains and 
    3-hydroxytetradecanoic acid units. They are a major constituent of the cell walls
    of Gram-negative bacteria.

Heuristic criteria used in this classifier:
  - At least one sugar-like ring is detected (rings of size 5–7 atoms that contain at least one oxygen and mostly carbons).
  - Presence of a carboxyl/carboxylic acid or ester motif is required. We look for a carbonyl bound to an oxygen ([CX3](=O)[OX2]).
  - Presence of a long aliphatic chain is expected. To capture this we look for at least eight contiguous carbon atoms 
    (as a proxy for at least 8 CH2 groups) using the SMARTS "CCCCCCCC".
  - Overall molecular weight is expected to be >1000 Da.

Note:
  This is a heuristic approach and some lipopolysaccharides may not meet all of these criteria perfectly.
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
    
    # Criterion 1: Check for sugar-like rings.
    # We count rings of size 5-7 that have at least one oxygen and nearly all atoms as carbons.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_ring_count = 0
    for ring in ring_info:
        if len(ring) not in [5, 6, 7]:
            continue
        # get element symbols for atoms in the ring
        symbols = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
        oxygen_count = symbols.count("O")
        carbon_count = symbols.count("C")
        # A typical sugar ring should have at least one oxygen and nearly all atoms are carbons.
        if oxygen_count >= 1 and carbon_count >= (len(ring) - 1):
            sugar_ring_count += 1
    # A relaxed threshold: require at least 1 such ring.
    if sugar_ring_count < 1:
        return False, f"Found only {sugar_ring_count} sugar-like ring(s); at least 1 expected"
    
    # Criterion 2: Check for a carboxyl/carboxylic acid or ester motif.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    if acid_pattern is None:
        return False, "Error in compiling acid SMARTS pattern"
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxyl/carboxylic acid motif ([CX3](=O)[OX2]) found"
    
    # Criterion 3: Look for a long aliphatic chain.
    # The previous SMARTS "[CH2]{8,}" failed to compile so we use "CCCCCCCC" as a heuristic to find eight contiguous carbons.
    fatty_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if fatty_chain_pattern is None:
        return False, "Error in compiling fatty chain SMARTS pattern"
    if not mol.HasSubstructMatch(fatty_chain_pattern):
        return False, "No long aliphatic chain (at least 8 contiguous carbons) found to represent fatty acid units"
    
    # Criterion 4: Check the molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a lipopolysaccharide"
    
    return True, ("Structure contains sugar-like ring(s), carboxyl/ester motifs, and long aliphatic chains "
                  "with a sufficiently high molecular weight, consistent with a lipopolysaccharide.")

# Example use (for testing purposes):
if __name__ == "__main__":
    test_smiles = "O=C(OC1C(O)C2(OC(C1OC3OC(C(O)C(C3O)O)COC(=O)/C=C/C=C/C=C/C)CO)OCC=4C2=C(O)C=C(O)C4)/C=C/C=C/CC(O)/C(=C/C=C/CCC(CC)C)/C"
    result, reason = is_lipopolysaccharide(test_smiles)
    print("Lipopolysaccharide classification:", result)
    print("Reason:", reason)