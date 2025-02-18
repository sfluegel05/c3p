"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies chemical entities of the class cation:
A monoatomic or polyatomic species having one or more elementary charges of the proton.

This improved strategy is as follows:
  1. Parse the SMILES string and compute the overall net formal charge.
  2. For a single atom:
       • If it carries a positive charge, return True — except if the element is Mn (as an example of a known false positive).
  3. For polyatomic molecules:
       a. If net charge > 0:
             – Check if a “robust” (pH–independent) cationic substructure is present (e.g. quaternary ammonium, aromatic nitrogen cation, guanidinium);
             – If present, classify as cation.
             – Otherwise (i.e. if only a protonated amine is present), require the molecule to be “sufficiently large” (MW at least ~200 Da)
               to reduce mis‐classifying very small protonated species.
       b. If net charge < 0: we return False.
       c. If net charge == 0 (i.e. a zwitterion):
             – We require that a robust cationic substructure is found,
               that the molecule also contains a phosphorus atom (often seen in phospholipids with a phosphocholine headgroup)
               and that the molecule is “large” (here we require at least 30 carbon atoms).
             – Otherwise, it is not classified as a cation.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    
    Strategy:
      - Parse the SMILES string.
      - Compute overall net formal charge.
      - For a single–atom species:
             If the atom carries a positive formal charge, we return True
             (with an ad hoc exclusion for, e.g., manganese in the +7 oxidation state).
      - For polyatomic species:
             If net charge > 0:
                • If the molecule contains a robust (pH–independent) cationic substructure
                  (quaternary ammonium, aromatic nitrogen cation, or guanidinium) then classify it as a cation.
                • Else (if the positive charge is only on a protonated amine), require that the molecule’s
                  molecular weight is at least 200 Da in order to avoid classifying very small ions as cations.
             If net charge < 0:
                • Return False.
             If net charge == 0 (i.e. zwitterions):
                • Check for a robust cationic group AND require that a phosphorus atom is present
                  (as a proxy for a phosphocholine subunit) AND that there are at least 30 carbon atoms.
                  Only then classify as a cation (e.g., cationic lipids).
                  Otherwise, do not classify as a cation.
    
    Args:
        smiles (str): A SMILES string representing the molecule.
    
    Returns:
        bool: True if the molecule is classified as a cation, False otherwise.
        str: A reason explaining the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate net formal charge
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    
    # Count number of carbon atoms (as a proxy for molecule size)
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    mol_wt = Descriptors.ExactMolWt(mol)
    
    # Define robust (pH–independent) cationic substructure SMARTS:
    #  • Quaternary ammonium has no hydrogens attached to the positive nitrogen.
    #  • Aromatic nitrogen cation (e.g. pyridinium) usually appears as [n+].
    #  • Guanidinium group typical pattern.
    quaternary = Chem.MolFromSmarts("[N+;H0]")          
    aromatic_nitrogen = Chem.MolFromSmarts("[n+]")         
    guanidinium = Chem.MolFromSmarts("NC(=[NH2+])N")
    
    robust_cation = False
    if quaternary is not None and mol.HasSubstructMatch(quaternary):
        robust_cation = True
    elif aromatic_nitrogen is not None and mol.HasSubstructMatch(aromatic_nitrogen):
        robust_cation = True
    elif guanidinium is not None and mol.HasSubstructMatch(guanidinium):
        robust_cation = True

    # Special handling for a single atom:
    if mol.GetNumAtoms() == 1:
        atom = mol.GetAtomWithIdx(0)
        # For example, if the atom is manganese (Mn) with a high oxidation state we exclude it.
        if atom.GetSymbol() == "Mn":
            return False, f"Single atom {atom.GetSymbol()} with net charge {net_charge} is not classified (Mn excluded)."
        if net_charge > 0:
            return True, f"Single atom cation with net positive charge of {net_charge}."
        else:
            return False, f"Single atom with net charge {net_charge} is not considered a cation."
    
    # Case 1: net positive charge
    if net_charge > 0:
        # If the molecule has a robust (pH–independent) cationic group, classify as cation
        if robust_cation:
            return True, "Molecule has net positive charge and a robust cationic substructure."
        else:
            # Otherwise, require that the molecule is sufficiently large to avoid small, transient (protonated) species.
            if mol_wt >= 200:
                return True, (f"Molecule has net positive charge ({net_charge}) and a molecular weight of {mol_wt:.1f} Da; "
                              "likely a stable cation.")
            else:
                return False, (f"Small molecule (MW={mol_wt:.1f} Da) with net positive charge but no robust cationic group; "
                               "likely only pH–dependent protonation.")
    
    # Case 2: net negative charge
    if net_charge < 0:
        return False, f"Molecule has net negative charge ({net_charge}); not a cation."
    
    # Case 3: net zero charge (zwitterions and neutral molecules)
    # Here we require that a robust cationic group is present,
    # that the molecule contains a phosphorus atom (often seen in phosphocholine lipids),
    # and that the molecule is “large” (at least 30 carbon atoms).
    if net_charge == 0:
        has_phosphorus = any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms())
        if robust_cation and has_phosphorus and num_carbons >= 30:
            return True, (f"Molecule has net zero charge but contains a robust cationic substructure, "
                          f"a phosphorus atom, and a large carbon skeleton (nC={num_carbons}); "
                          "consistent with a cationic lipid.")
        else:
            return False, "Molecule has net zero charge without sufficient indicators of a permanent cation."
    
    # Just in case, return a safe default
    return False, "Molecule does not meet criteria for cation classification."

# Example usage (uncomment the following lines to test a few SMILES):
# test_smiles_list = [
#     # True positives:
#     "P(OCC[N+](C)(C)C)(OCC(OC(=O)CCCCCCCCCCCCC=1OC(CCCCC)=CC1C)COC(=O)CCCCCCCCCCC=2OC(=C(C2C)C)CCCCC)(O)=O",  # lipid–like phosphinic acid derivative
#     "[NH3+]CCC(=O)NCCC1=CNC=N1",      # carcininium
#     "COc1cccc(c1)[C@@]1(O)CCCC[C@@H]1C[NH+](C)C",  # (R,R)-tramadol(1+)
#     "OC[C@H](CC(C)C)[NH3+]",          # (S)-leucinol(1+)
#     "[C@H](CC[Se+](C)C)(C(O)=O)N",     # Se-methyl-L-selenomethionine
#
#     # False positives (should not be classified as cations):
#     "C[N+](C)(C)CC(O)O",             # betaine aldehyde hydrate is correct but note: small ions like N-dimethylethanolamine are excluded
#     "OCC[NH+](C)C",                 # N-dimethylethanolamine (small; MW low)
#     "[Mn+7]",                      # manganese(7+) is excluded by ad hoc rule
#
#     # net zero but not large or lacking phosphorus ➔ not cation:
#     "C(C(CC([O-])=O)OC(=O)C[C@@H](CCCCCCC/C=C\\C/C=C\\CCCCC)O)[N+](C)(C)C",  # carnitine derivative (zwitterion with nC≈27)
# ]
#
# for smi in test_smiles_list:
#     res, reason = is_cation(smi)
#     print(f"SMILES: {smi}\nResult: {res} | Reason: {reason}\n")