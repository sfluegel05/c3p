"""
Classifies: CHEBI:32878 alkene
"""
"""
Classifies: Alkene (acyclic hydrocarbon having exactly one C=C double bond
with the general formula CnH2n)
"""
import re
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene (acyclic mono‐unsaturated hydrocarbon)
    based on its SMILES string. The definition used is that the molecule must:
      1. Be a valid molecule with only carbon and hydrogen.
      2. Be acyclic (no rings).
      3. Have exactly one carbon–carbon double bond.
      4. Have a molecular formula that exactly fits CnH2n.
      
    Note: If the molecule does not meet these strict criteria then it is rejected.
    (For ambiguous cases additional heuristics may be needed. If the task appears too hard,
    one may return (None, None).)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        tuple: (bool, str) where the boolean indicates whether the molecule qualifies as an alkene,
               and the string gives the reason.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that the molecule is strictly a hydrocarbon (only C and H)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6):
            return False, "Molecule contains atoms other than carbon and hydrogen"

    # Check that the molecule is acyclic (has no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic; expected an acyclic structure"

    # Use a SMARTS pattern that explicitly looks for a C=C bond where both carbons are acyclic.
    # (The pattern [C;!R]=[C;!R] ensures both carbons are not in a ring.)
    alkene_smarts = "[C;!R]=[C;!R]"
    alkene_query = Chem.MolFromSmarts(alkene_smarts)
    double_bond_matches = mol.GetSubstructMatches(alkene_query)
    dbl_bond_count = len(double_bond_matches)
    if dbl_bond_count != 1:
        return False, f"Molecule has {dbl_bond_count} C=C bond(s); expected exactly 1"

    # Compute the molecular formula using RDKit’s built‐in function.
    formula = rdMolDescriptors.CalcMolFormula(mol)
    # Use a regex to extract the carbon and hydrogen counts.
    match = re.match(r"^C(\d*)H(\d+)$", formula)
    if not match:
        return False, f"Formula {formula} did not match expected pattern CnH2n"
    # If the number is missing after C assume 1 (i.e. "CH4" -> 1 carbon)
    c_count = int(match.group(1)) if match.group(1) != "" else 1
    h_count = int(match.group(2))
    
    # Check that the formula exactly follows CnH2n
    if h_count != 2 * c_count:
        return False, f"Molecular formula is {formula}, which does not match CnH2n"
        
    # Passed all tests: classify as a valid alkene.
    return True, f"Molecule is a valid acyclic alkene with 1 C=C bond and formula {formula}"

# If the task proves too ambiguous based on further tests one might instead fallback:
# return None, None