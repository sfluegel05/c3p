"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: Alkene
Definition: An acyclic branched or unbranched hydrocarbon having exactly one carbon–carbon double bond 
with the general formula CnH2n (mono‐unsaturated).
Any molecule that fails one of the following criteria is rejected:
  1. A valid molecule.
  2. Contains only C and H atoms.
  3. Is acyclic (has no rings).
  4. Contains exactly one C=C double bond (between two carbon atoms).
  5. Has a molecular formula that exactly matches CnH2n.
  
Improvement: The molecule is first sanitized and explicit hydrogens are removed so that the computed molecular 
formula is more reliable even when stereochemistry or explicit hydrogens are provided.
"""
import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene (acyclic mono‐unsaturated hydrocarbon)
    based on its SMILES string. Specifically, the molecule must:
      1. Be a valid molecule.
      2. Contain only carbon and hydrogen atoms.
      3. Be acyclic.
      4. Have exactly one C=C double bond (between carbons).
      5. Have a molecular formula of the form CnH2n.
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) where the boolean is True if the molecule qualifies as an alkene, and
               the string provides a reason for the classification outcome.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize and remove explicit hydrogens so that the molecular formula is computed uniformly.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"
    mol = Chem.RemoveHs(mol)
    
    # Check that the molecule contains only carbon (atomic number 6) and hydrogen (atomic number 1)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):
            return False, "Molecule contains atoms other than carbon and hydrogen"
    
    # Verify the molecule is acyclic (i.e. it has no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic; expected an acyclic structure"
    
    # Count the number of C=C (carbon–carbon double bonds).
    cc_double_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                cc_double_count += 1
    if cc_double_count != 1:
        return False, f"Molecule has {cc_double_count} C=C bond(s); expected exactly 1"
    
    # Compute the molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    
    # Match the formula with a strict pattern “CnHn” where hydrogen count must equal 2 * carbon count.
    # Using fullmatch ensures the formula contains only C and H in that order.
    match = re.fullmatch(r"C(\d*)H(\d+)", formula)
    if not match:
        return False, f"Formula {formula} did not match expected pattern CnH2n"
    
    # If the carbon count group is empty, that means there is one carbon.
    c_count = int(match.group(1)) if match.group(1) != "" else 1
    h_count = int(match.group(2))
    
    if h_count != 2 * c_count:
        return False, f"Molecular formula is {formula}, which does not match CnH2n"
    
    return True, f"Molecule is a valid acyclic alkene with 1 C=C bond and formula {formula}"