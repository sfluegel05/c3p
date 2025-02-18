"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: Volatile Organic Compound (VOC)
Definition: “Any organic compound having an initial boiling point less than or equal to 250 °C 
(measured at the standard atmospheric pressure of 101.3 kPa).”
Heuristic strategy:

  1. The molecule must contain carbon.
  2. Compute molecular descriptors: molecular weight (MW) and topological polar surface area (TPSA).
  3. Simple molecules (defined as having no rings and at most one non‐halogen heteroatom)
     are allowed a higher MW cutoff (350 Da) because many unfunctionalized long‐chain alcohols are VOC even if MW >300.
     For “more complex” molecules, we require MW <= 300.
  4. If the molecule contains too many heteroatoms (i.e. >2 non‐halogen heteroatoms) we assume the molecule is too functionalized.
  5. Additionally, if the molecule shows an ester group (and not a carboxylic acid) we assume it is not a VOC,
     as many esters have higher boiling points despite low MW/TPSA.
  6. Finally, if MW is less than the cutoff and TPSA is less than 60 Å² then we classify it as VOC.
     
Note: This is a rough heuristic and does not constitute a full QSPR prediction.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on a heuristic estimation.
    
    For this heuristic:
      - All molecules must be organic (contain at least one carbon atom).
      - We calculate molecular weight (MW) and topological polar surface area (TPSA).
      - For relatively simple molecules (no rings and at most 1 non‐halogen heteroatom), the MW cutoff is raised to 350 Da.
        For more complex molecules the cutoff remains 300 Da.
      - Molecules with >2 non‐halogen heteroatoms are assumed to be too highly functionalized.
      - In addition, if an ester group is present (and no acid group is present) the compound is not considered VOC.
      - Finally, we require TPSA to be below 60 Å².
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the compound is classified as VOC, False otherwise.
        str: Explanation string.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic check: the molecule must be organic (contain at least one carbon).
    if not any(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms()):
        return False, "Not organic (contains no carbon)"
    
    # Calculate molecular weight and TPSA.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    
    # Count non-halogen heteroatoms (exclude C, H, and common halogens: F, Cl, Br, I).
    allowed_halogens = {9, 17, 35, 53}
    hetero_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6) and atom.GetAtomicNum() not in allowed_halogens)
    
    # Check if the molecule is "simple": no rings and at most 1 non-halogen heteroatom.
    simple = (mol.GetRingInfo().NumRings() == 0) and (hetero_count <= 1)
    
    # Set molecular weight cutoff.
    cutoff_mw = 350 if simple else 300
    
    # Check for excessive functionalization: if too many heteroatoms, reject.
    if hetero_count > 2:
        return False, (f"Too many heteroatoms ({hetero_count} non-halogen atoms), suggesting high functionality and likely a high boiling point.")
    
    # Use SMARTS to detect ester groups: [CX3](=O)[OX2H0] (ester) but not carboxylic acids [CX3](=O)[OX1H]
    ester_smarts = "[CX3](=O)[OX2H0]"
    acid_smarts = "[CX3](=O)[OX1H]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    has_acid = mol.HasSubstructMatch(acid_pattern)
    if has_ester and not has_acid:
        return False, "Contains an ester group (without acid functionality), which tends to increase the boiling point."
    
    # Final heuristic: check if MW and TPSA meet criteria.
    if mol_wt <= cutoff_mw and tpsa < 60:
        return True, (f"Estimated as VOC: MW ({mol_wt:.1f} Da) <= {cutoff_mw} and TPSA ({tpsa:.1f} Å²) <60, "
                      "suggesting a low boiling point (<=250 °C).")
    else:
        return False, (f"Estimated not VOC: MW ({mol_wt:.1f} Da) and TPSA ({tpsa:.1f} Å²) do not meet criteria "
                       f"for low boiling point (<=250 °C) with cutoff MW = {cutoff_mw} Da.")

# Example usage:
if __name__ == "__main__":
    # A few test SMILES examples including both true and false positives/negatives
    test_cases = {
        "nonan-2-ol": "CCCCCCCC(C)O",
        "decan-2-ol": "CCCCCCCCC(C)O",
        "2-dodecene": "[H]C(C)=C([H])CCCCCCCCC",
        "henicosan-3-ol": "CCCCCCCCCCCCCCCCCC(O)CC",      # previously false negative: MW ~312
        "hexacosan-4-ol": "CCCCCCCCCCCCCCCCCCCCCCCC(O)CCC", # previously false negative: MW ~382
        "pramocaine": "CCCCOC1=CC=C(C=C1OCCCN2CCOCC2)C",
        "thiophene": "c1ccsc1"
    }
    
    for name, sm in test_cases.items():
        is_voc, reason = is_volatile_organic_compound(sm)
        print(f"Name: {name}\n  SMILES: {sm}\n  VOC: {is_voc}\n  Reason: {reason}\n")