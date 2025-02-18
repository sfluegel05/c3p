"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies chemical entities of the class cation:
A monoatomic or polyatomic species having one or more elementary charges of the proton.

Our improved strategy:
  1. Parse the SMILES and compute the overall net formal charge.
  2. If the molecule is a single atom and its formal charge > 0, return True.
  3. If net charge > 0, classify as a cation.
  4. If net charge < 0, not a cation.
  5. For net zero molecules (i.e. potential zwitterions) we look for a robust, pH‐independent cationic
     substructure (quaternary ammonium, aromatic nitrogen or guanidinium) AND check if the molecule
     is “large” (has at least 25 carbon atoms). This extra size filter is intended to flag lipids
     (such as phosphatidylcholines) that, although formally neutral, behave as cationic species.
     
Any molecule that does not meet these criteria is not classified as a cation.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    
    Strategy:
      - Parse the SMILES and compute overall net formal charge.
      - If the molecule is a single atom with positive charge, classify as a cation.
      - If net formal charge is >0 then classify as cation.
      - If net formal charge is <0 then not a cation.
      - If overall net charge is 0, check for the presence of a robust (pH‐independent) cationic substructure:
           * For our purposes we define these as:
                 - quaternary ammonium: SMARTS "[N+;H0]"
                 - aromatic nitrogen cation (e.g. pyridinium): "[n+]"
                 - guanidinium group: "NC(=[NH2+])N"
        If one of these is found AND the molecule is “large” (we require at least 25 carbon atoms)
        then we classify it as a cation (e.g. many phosphatidylcholine lipids). Otherwise, we do not.
    
    Args:
        smiles (str): A SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule is classified as a cation, False otherwise.
        str: A reason explaining the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute overall net formal charge (sum of formal charges on all atoms)
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    
    # Special case: if molecule is a single atom, treat positive atoms as cations.
    if mol.GetNumAtoms() == 1:
        if net_charge > 0:
            return True, f"Single atom cation with net positive charge of {net_charge}."
        else:
            return False, f"Single atom with net charge {net_charge} is not considered a cation."
    
    mol_wt = Descriptors.ExactMolWt(mol)
    # Count number of carbon atoms as a quick estimate of molecular size.
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Define SMARTS patterns for robust, pH‐independent cationic groups.
    quaternary = Chem.MolFromSmarts("[N+;H0]")          # e.g., [N+](C)(C)C (no hydrogen attached)
    aromatic_nitrogen = Chem.MolFromSmarts("[n+]")         # e.g., pyridinium
    guanidinium = Chem.MolFromSmarts("NC(=[NH2+])N")        # typical guanidinium group
    
    # ----------- Case 1: net_charge > 0 -----------
    if net_charge > 0:
        return True, f"Molecule has net positive charge ({net_charge}); classified as cation."
    
    # ----------- Case 2: net_charge < 0 -----------
    if net_charge < 0:
        return False, f"Molecule has net negative charge ({net_charge}); not a cation."
    
    # ----------- Case 3: net_charge == 0 -----------
    # At net zero overall, many molecules are zwitterionic.
    # We then check if the molecule contains a robust, permanent cationic substructure.
    has_robust_cation = False
    if quaternary is not None and mol.HasSubstructMatch(quaternary):
        has_robust_cation = True
    elif aromatic_nitrogen is not None and mol.HasSubstructMatch(aromatic_nitrogen):
        has_robust_cation = True
    elif guanidinium is not None and mol.HasSubstructMatch(guanidinium):
        has_robust_cation = True
    
    if has_robust_cation:
        # Use number of carbons as a proxy for a large scaffold (e.g., lipid chains).
        if num_carbons >= 25:
            return True, (f"Molecule has net zero charge but contains a robust "
                          f"cationic functional group and a large carbon skeleton "
                          f"(nC={num_carbons}), consistent with a cationic lipid.")
        else:
            return False, (f"Molecule has net zero charge and a robust cationic functional "
                           f"group but is small (nC={num_carbons}); likely its positive site is "
                           "only pH-dependent and is balanced by an anionic group.")
    
    return False, f"Molecule has net zero charge with no unambiguous permanent cationic substructure."

# Example usage (uncomment to test):
# test_smiles_list = [
#     # True positives:
#     "P(OCC[N+](C)(C)C)(OCC(OC(=O)CCCCCCCCCCCCC=1OC(CCCCC)=CC1C)COC(=O)CCCCCCCCCCC=2OC(=C(C2C)C)CCCCC)(O)=O",  # lipid-like phosphinic acid derivative
#     "[NH3+]CCC(=O)NCCC1=CNC=N1",  # carcininium with +1 charge
#     "COc1cccc(c1)[C@@]1(O)CCCC[C@@H]1C[NH+](C)C",  # (R,R)-tramadol(1+)
#     "OC[C@H](CC(C)C)[NH3+]",  # (S)-leucinol(1+)
#     "[C@H](CC[Se+](C)C)(C(O)=O)N",  # Se-methyl-L-selenomethionine
#
#     # False positives (should be classified as not cations):
#     "C(C(CC([O-])=O)OC(=O)C[C@@H](CCCCCCC/C=C\\C/C=C\\CCCCC)O)[N+](C)(C)C",  # a carnitine derivative (zwitterion)
#     "[Mn+7]",  # manganese +7 (exotic, but net +7 on an atom – note: our rule would classify single atom positive, so such cases remain ions)
# ]
#
# for smi in test_smiles_list:
#     res, reason = is_cation(smi)
#     print(f"SMILES: {smi}\nResult: {res} | Reason: {reason}\n")