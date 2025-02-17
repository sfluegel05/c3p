"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: Volatile Organic Compound (VOC)
Definition: “Any organic compound having an initial boiling point less than or equal 
to 250 °C (482 °F) measured at a standard atmospheric pressure of 101.3 kPa.”

Revised Heuristic Strategy:
  1. The molecule must contain at least one carbon atom.
  2. Calculate descriptors: molecular weight (MW) and topological polar surface area (TPSA).
  3. Count non‐halogen heteroatoms (excluding H, carbon and common halogens F, Cl, Br, I).
  4. Assess ring count – if acyclic, allow higher MW cutoff (<=400 Da); if one or more rings,
     use a stricter cutoff (<=300 Da).
  5. If the molecule is “highly functionalized” (≥3 non‐halogen heteroatoms, or if it has 
     ≥2 rings and ≥2 heteroatoms) then assume it has a high boiling point.
  6. If the molecule contains a problematic functional group (an ester group without also 
     having a carboxylic acid group, a carboxylic acid, aromatic hydroxyl (phenol), aromatic 
     amine, or a cyclic ketone), disqualify as VOC.
  7. Finally, if MW is below the cutoff and TPSA is <60 Å², classify as VOC.

Note: This heuristic is approximate.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on a 
    heuristic estimation.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as VOC, False otherwise.
        str: Explanation of the classification.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is organic (has at least one carbon atom)
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not organic (contains no carbon)"
    
    # Calculate descriptors: molecular weight (MW) and topological polar surface area (TPSA)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    
    # Count the non-halogen heteroatoms (exclude H=1, C=6 and common halogens: F=9, Cl=17, Br=35, I=53).
    allowed_halogens = {9, 17, 35, 53}
    hetero_count = sum(1 for atom in mol.GetAtoms() 
                        if atom.GetAtomicNum() not in (1, 6) and atom.GetAtomicNum() not in allowed_halogens)
    
    # Get ring information: number of rings
    ring_count = mol.GetRingInfo().NumRings()
    # Define "simple" as having no rings.
    simple = (ring_count == 0)
    
    # Set molecular weight cutoff based on ring content:
    # For simple (acyclic) molecules, allow a higher cutoff (<=400 Da);
    # for molecules with rings, use a stricter cutoff (<=300 Da).
    cutoff_mw = 400 if simple else 300

    # Extra rejection based on high functionality:
    if hetero_count >= 3:
        return False, f"Too many heteroatoms ({hetero_count}), suggesting high functionality and high boiling point."
    if (not simple) and hetero_count >= 2 and ring_count >= 2:
        return False, f"Multiple rings ({ring_count}) and heteroatoms ({hetero_count}) suggest non‐volatile functionality."
    
    # Define SMARTS patterns for problematic functional groups:
    # 1. Ester group (but allow if there is also a carboxylic acid).
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0]")
    acid_pattern  = Chem.MolFromSmarts("[CX3](=O)[OX1H]")
    if mol.HasSubstructMatch(ester_pattern) and not mol.HasSubstructMatch(acid_pattern):
        return False, "Contains an ester group (without accompanying acid functionality), which tends to increase the boiling point."
    
    # 2. Carboxylic acid group.
    if mol.HasSubstructMatch(acid_pattern):
        return False, "Contains a carboxylic acid group, which tends to raise the boiling point."
    
    # 3. Aromatic hydroxyl (phenol): oxygen directly attached to an aromatic carbon.
    phenol_pattern = Chem.MolFromSmarts("c[OH]")
    if mol.HasSubstructMatch(phenol_pattern):
        return False, "Contains an aromatic hydroxyl group (phenol), known to increase boiling point."
    
    # 4. Aromatic primary amine: aromatic carbon with attached NH2.
    arylamine_pattern = Chem.MolFromSmarts("c[NH2]")
    if mol.HasSubstructMatch(arylamine_pattern):
        return False, "Contains an aromatic amine group, which tends to increase boiling point."
    
    # 5. Cyclic ketone: carbonyl group within a ring.
    cyclic_ketone = Chem.MolFromSmarts("[R][CX3](=O)[R]")
    if mol.HasSubstructMatch(cyclic_ketone):
        return False, "Contains a cyclic ketone group, which tends to raise the boiling point."
        
    # Final decision based on descriptors:
    # We require TPSA to be under 60 Å² and molecular weight to be within the cutoff.
    if mol_wt <= cutoff_mw and tpsa < 60:
        return True, (f"Estimated as VOC: MW ({mol_wt:.1f} Da) <= {cutoff_mw} and TPSA ({tpsa:.1f} Å²) < 60, "
                      "suggesting a low boiling point (<=250 °C).")
    else:
        return False, (f"Estimated not VOC: MW ({mol_wt:.1f} Da) and TPSA ({tpsa:.1f} Å²) do not meet criteria "
                       f"for low boiling point (<=250 °C) with cutoff MW = {cutoff_mw} Da.")

# Example usage (for testing)
if __name__ == "__main__":
    test_cases = {
        "nonan-2-ol": "CCCCCCCC(C)O",
        "decan-2-ol": "CCCCCCCCC(C)O",
        "2-dodecene": "[H]C(C)=C([H])CCCCCCCCC",
        "henicosan-3-ol": "CCCCCCCCCCCCCCCCCC(O)CC",
        "hexacosan-4-ol": "CCCCCCCCCCCCCCCCCCCCCCCC(O)CCC",
        "thiophene": "c1ccsc1",
        "4-methylpent-3-enoic acid": "OC(=O)CC=C(C)C",
        "pre-flavunoidine(2+)": "C1(C2CC34C(C1(OC2(CC3)C)C)CCC4C)[NH2+]CCCCC[NH+](C)C",
    }
    
    for name, sm in test_cases.items():
        voc, reason = is_volatile_organic_compound(sm)
        print(f"Name: {name}\n  SMILES: {sm}\n  VOC: {voc}\n  Reason: {reason}\n")