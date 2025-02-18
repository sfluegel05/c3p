"""
Classifies: CHEBI:134179 volatile organic compound
"""
"""
Classifies: Volatile Organic Compound (VOC)
Definition: “Any organic compound having an initial boiling point less than or equal to 250 °C 
(measured at a standard atmospheric pressure of 101.3 kPa).”
This implementation uses a heuristic method as a surrogate for actual boiling point prediction.
Heuristic:
  - The molecule must be organic (contain carbon).
  - We calculate the molecular weight (MW) and topological polar surface area (TPSA).
  - Molecules with MW <= 300 Daltons and TPSA < 60 Å² are assumed to have an initial boiling
    point likely <= 250 °C.
Note: This is only a rough estimate and not a substitute for a true QSPR model.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_volatile_organic_compound(smiles: str):
    """
    Determines if a molecule is a volatile organic compound (VOC) based on a heuristic estimation.
    For our purposes, we assume that an organic molecule with a molecular weight <= 300 Da and a 
    topological polar surface area < 60 Å² will have an initial boiling point <= 250 °C.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the compound is classified as a VOC, False otherwise.
        str: Explanation of the reasoning.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic check: ensure the molecule is organic (contains at least one carbon)
    if not any(atom.GetSymbol() == "C" for atom in mol.GetAtoms()):
        return False, "Not an organic compound (contains no carbon)"

    # Calculate molecular descriptors.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)

    # Heuristic interpretation:
    # Many VOCs are relatively small and nonpolar. Here we classify molecules with MW <=300 Da and
    # TPSA < 60 Å² as likely to have a boiling point <= 250 °C.
    if mol_wt <= 300 and tpsa < 60:
        return True, (f"Estimated as VOC: MW ({mol_wt:.1f} Da) <=300 and TPSA ({tpsa:.1f} Å²) <60, "
                      "suggesting a low boiling point (<=250 °C).")
    else:
        return False, (f"Estimated not VOC: MW ({mol_wt:.1f} Da) and TPSA ({tpsa:.1f} Å²) do not meet criteria "
                       "for low boiling point (<=250 °C).")
                        
# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "CCCCCCCC(C)O",       # nonan-2-ol
        "CCCCCCCCC(C)O",       # decan-2-ol
        "[H]C(C)=C([H])CCCCCCCCC",  # 2-dodecene
        "c1ccsc1",            # thiophene
        "CCCOC(=O)C1CCCCC1"    # ethyl cyclohexanecarboxylate (example)
    ]
    for sm in test_smiles:
        classification, reason = is_volatile_organic_compound(sm)
        print(f"SMILES: {sm}\n  VOC: {classification}\n  Reason: {reason}\n")