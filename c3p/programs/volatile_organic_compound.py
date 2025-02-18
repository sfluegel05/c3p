"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: Volatile Organic Compound (VOC)
Definition: “Any organic compound having an initial boiling point less than or equal 
to 250 °C (482 °F) measured at a standard atmospheric pressure of 101.3 kPa.”
Revised Heuristic Strategy:
  1. The molecule must contain carbon.
  2. Calculate descriptors: molecular weight (MW) and topological polar surface area (TPSA).
  3. Count non‐halogen heteroatoms.
  4. If the molecule is “simple” (no rings), allow a higher MW cutoff (350 Da). 
     Otherwise (one or more rings) use a stricter cutoff (300 Da).
  5. If the molecule is too functionalized (≥3 non‐halogen heteroatoms, or is “complex”
     with 2 or more heteroatoms and ≥2 rings) then assume it has a high boiling point.
  6. Additionally, if it contains an ester group (without an accompanying acid group),
     a carboxylic acid group, an aromatic hydroxyl (phenol) or aromatic amine, or 
     a cyclic ketone substructure, then it is assumed not to be VOC.
  7. Finally, if MW and TPSA are under the cutoff values (and TPSA < 60 Å²) classify as VOC.
     
Note: This heuristic is still approximate.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.Crippen import MolLogP

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on a 
    heuristic estimation.
    
    The heuristic:
      - Requires the molecule be organic (contain at least one carbon atom).
      - Calculates molecular weight (MW) and topological polar surface area (TPSA).
      - Considers “simple” molecules (with no rings) to tolerate higher MW (<=350 Da),
        while “complex” molecules (with rings) must have MW <=300 Da.
      - Additionally, if the molecule is highly functionalized (≥3 non‐halogen heteroatoms
        or (≥2 heteroatoms and ≥2 rings)) it is assumed to have a high boiling point.
      - Further, if the molecule shows problematic groups (ester without acid, a carboxylic acid,
        aromatic hydroxyl or arylamine substituents, or a cyclic ketone group) it is disqualified.
      - Finally, if TPSA is under 60 Å² and the descriptor values are below cutoff, it is classified as VOC.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as VOC, False otherwise.
        str: Explanation string.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ensure the molecule is organic (contains at least one carbon)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not organic (contains no carbon)"
    
    # Compute descriptors
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    
    # Count non-halogen heteroatoms (exclude C, H and common halogens: F, Cl, Br, I)
    allowed_halogens = {9, 17, 35, 53}
    hetero_count = sum(1 for atom in mol.GetAtoms() 
                       if atom.GetAtomicNum() not in (1,6) and atom.GetAtomicNum() not in allowed_halogens)
    
    # Get ring information: count how many rings are present.
    ring_count = mol.GetRingInfo().NumRings()
    
    # Determine if molecule is "simple" (aliphatic, no rings)
    simple = (ring_count == 0)
    
    # Set molecular weight cutoff: simple: <=350 Da; complex: <=300 Da.
    cutoff_mw = 350 if simple else 300
    
    # Extra rejection if too many heteroatoms: if ≥3 non‐halogen heteroatoms
    if hetero_count >= 3:
        return False, f"Too many heteroatoms ({hetero_count}), suggesting high functionality and high boiling point."
    # Also if molecule is not simple and has 2 or more heteroatoms with ≥2 rings, reject.
    if (not simple) and hetero_count >=2 and ring_count >=2:
        return False, f"Multiple rings and heteroatoms ({ring_count} rings, {hetero_count} heteroatoms) suggest non‐volatile functionality."
    
    # SMARTS patterns for problematic functional groups:
    # Ester group (but allow if also a carboxylic acid) -> non VOC.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    acid_pattern  = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    if mol.HasSubstructMatch(ester_pattern) and not mol.HasSubstructMatch(acid_pattern):
        return False, "Contains an ester group (without acid functionality), which tends to increase the boiling point."
    
    # Carboxylic acid (even if present with ester, acids tend to raise boiling point)
    if mol.HasSubstructMatch(acid_pattern):
        return False, "Contains a carboxylic acid group, which tends to raise the boiling point."
    
    # Aromatic hydroxyl (phenol) – detect oxygen bound directly to an aromatic carbon.
    phenol_pattern = Chem.MolFromSmarts("c[OH]")
    if mol.HasSubstructMatch(phenol_pattern):
        return False, "Contains an aromatic hydroxyl group (phenol), known to increase boiling point."
    
    # Aromatic primary amine – detect pattern: aromatic carbon attached to NH2.
    arylamine_pattern = Chem.MolFromSmarts("c[NH2]")
    if mol.HasSubstructMatch(arylamine_pattern):
        return False, "Contains an aromatic amine group, which tends to increase the boiling point."
    
    # Cyclic ketone: carbonyl group within a ring
    cyclic_ketone = Chem.MolFromSmarts("[R][CX3](=O)[R]")
    if mol.HasSubstructMatch(cyclic_ketone):
        return False, "Contains a cyclic ketone group, which tends to raise the boiling point."
    
    # Final decision based solely on descriptors.
    # Also require TPSA to be less than 60 Å².
    if mol_wt <= cutoff_mw and tpsa < 60:
        return True, (f"Estimated as VOC: MW ({mol_wt:.1f} Da) <= {cutoff_mw} and TPSA ({tpsa:.1f} Å²) < 60, "
                      "suggesting a low boiling point (<=250 °C).")
    else:
        return False, (f"Estimated not VOC: MW ({mol_wt:.1f} Da) and TPSA ({tpsa:.1f} Å²) do not meet criteria "
                       f"for low boiling point (<=250 °C) with cutoff MW = {cutoff_mw} Da.")

# Example usage (for manual testing)
if __name__ == "__main__":
    test_cases = {
        "nonan-2-ol": "CCCCCCCC(C)O",
        "decan-2-ol": "CCCCCCCCC(C)O",
        "2-dodecene": "[H]C(C)=C([H])CCCCCCCCC",
        "henicosan-3-ol": "CCCCCCCCCCCCCCCCCC(O)CC",
        "hexacosan-4-ol": "CCCCCCCCCCCCCCCCCCCCCCCC(O)CCC",
        "thiophene": "c1ccsc1",
        "4-(2-aminopropyl)phenol": "CC(CC1=CC=C(C=C1O)N)O",  # expected to be rejected
        "piperitone": "CC(C)C1CCC(C)=CC1=O",                # cyclic ketone; reject
    }
    
    for name, sm in test_cases.items():
        voc, reason = is_volatile_organic_compound(sm)
        print(f"Name: {name}\n  SMILES: {sm}\n  VOC: {voc}\n  Reason: {reason}\n")