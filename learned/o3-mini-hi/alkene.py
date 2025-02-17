"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: Alkene (acyclic hydrocarbon having exactly one C=C bond with the general formula CnH2n)
An alkene here must:
  1. Be a valid molecule.
  2. Contain only carbon and hydrogen.
  3. Be acyclic (contain no rings).
  4. Contain exactly one carbon–carbon double bond.
  5. Have a molecular formula that exactly matches CnH2n.
Any molecule failing one of these criteria is rejected.
"""
import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene (acyclic mono‐unsaturated hydrocarbon)
    based on its SMILES string. Specifically, the molecule must:
      1. contain only carbon and hydrogen (no heteroatoms),
      2. be acyclic,
      3. have exactly one C=C (carbon–carbon double bond) and
      4. have a molecular formula of the form CnH2n.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        tuple: (bool, str) where the boolean is True if the molecule qualifies as an alkene and
               the string gives the reason; otherwise False with an explanation.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is strictly a hydrocarbon (only C and H atoms)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):
            return False, "Molecule contains atoms other than carbon and hydrogen"
    
    # Check that the molecule is acyclic (has no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic; expected an acyclic structure"
    
    # Count the number of carbon–carbon double bonds.
    # We iterate over each bond once and count those bonds that are:
    #  - a double bond, and
    #  - between two carbon atoms.
    cc_double_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                cc_double_count += 1
    
    if cc_double_count != 1:
        return False, f"Molecule has {cc_double_count} C=C bond(s); expected exactly 1"
    
    # Compute the molecular formula using RDKit
    formula = rdMolDescriptors.CalcMolFormula(mol)
    # Use a regex to extract the number of carbons and hydrogens.
    # This expects the formula in the order C...H... (e.g., C11H22)
    match = re.match(r"^C(\d*)H(\d+)$", formula)
    if not match:
        return False, f"Formula {formula} did not match expected pattern CnH2n"
    
    # If the carbon count is missing (as in "CH4") assume one carbon.
    c_count = int(match.group(1)) if match.group(1) != "" else 1
    h_count = int(match.group(2))
    
    # Check that the formula exactly follows CnH2n for an acyclic mono‐unsaturated hydrocarbon.
    if h_count != 2 * c_count:
        return False, f"Molecular formula is {formula}, which does not match CnH2n"
    
    return True, f"Molecule is a valid acyclic alkene with 1 C=C bond and formula {formula}"