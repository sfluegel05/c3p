"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies chemical entities of the class cation:
A monoatomic or polyatomic species having one or more elementary charges of the proton.
This improved classifier combines overall net charge, molecular weight and SMARTS matching for
permanent cationic functional groups to better separate standalone cations from zwitterions
or molecules with only pH‐dependent protonations.
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    
    Improved strategy:
      1. Compute the overall net formal charge.
      2. For net positive species:
           a. If the molecule is a single atom (e.g. metal cation) => classify as cation.
           b. If the molecular weight is low (<150 Da), require the presence of a robust 
              (permanent) cationic substructure (e.g. quaternary ammonium, aromatic nitrogen, 
              guanidinium).
           c. For heavier molecules (>=150 Da), also reject those that combine a protonated amine 
              with a carboxylate (suggesting a zwitterion) unless a permanent cationic group is present.
      3. For net zero species, check if the sum of positive charges outweighs the negatives 
         or if a permanent cationic fragment is present.
      4. Net negative species are not cations.
      
    Returns:
        bool: True if classified as a cation, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute overall net formal charge
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    mol_wt = Descriptors.ExactMolWt(mol)
    
    # Define SMARTS patterns for permanent cationic groups:
    # Quaternary ammonium (no attached hydrogen – pH‐independent)
    quaternary = Chem.MolFromSmarts("[N+;H0]")
    # Aromatic nitrogen cation (e.g. pyridinium)
    aromatic_nitrogen = Chem.MolFromSmarts("[n+]")
    # Guanidinium group (commonly seen in arginine derivatives)
    guanidinium = Chem.MolFromSmarts("NC(=[NH2+])N")
    # Generic positive nitrogen (any nitrogen with positive formal charge)
    positive_n = Chem.MolFromSmarts("[N+]")
    # (For our purposes we treat a positively charged N that is not quaternary/aromatic/guanidinium
    #  as pH-dependent and subject to molecular-weight filtering.)
    
    # Define a pattern for a carboxylate group (which may indicate partners in zwitterions)
    carboxylate = Chem.MolFromSmarts("C(=O)[O-]")
    
    # Special check: if the molecule is a single atom and has a net positive charge,
    # classify it as a cation (e.g. [In+], [Fe++], [Cr+6])
    if mol.GetNumAtoms() == 1:
        if net_charge > 0:
            return True, f"Single atom cation with net positive charge of {net_charge}."
        else:
            return False, f"Single atom with net charge {net_charge} is not a cation."
    
    # ---------------------------
    # Case 1. net_charge > 0
    if net_charge > 0:
        # Helper: does the molecule have any robust (pH-independent) cationic substructure?
        has_permanent = (mol.HasSubstructMatch(quaternary) or 
                         mol.HasSubstructMatch(aromatic_nitrogen) or 
                         (guanidinium is not None and mol.HasSubstructMatch(guanidinium)))
        
        # For smaller molecules, only count permanent groups.
        if mol_wt < 150:
            if has_permanent:
                return True, f"Small molecule (MW={mol_wt:.1f} Da) has net positive charge {net_charge} and contains a permanent cationic group."
            else:
                return False, f"Small molecule (MW={mol_wt:.1f} Da) has net positive charge {net_charge} but lacks a robust cationic substructure (likely only pH-dependent protonation)."
        else:
            # For heavier species, also allow a general positive site—but be careful about zwitterions.
            # If the only positive site is pH-dependent (i.e. positive nitrogen not matching a permanent group)
            # and there is a carboxylate present, we suspect a zwitterion so do not classify as a cation.
            has_positive = mol.HasSubstructMatch(positive_n)
            has_carboxylate = mol.HasSubstructMatch(carboxylate)
            if not has_permanent:
                if has_positive and has_carboxylate:
                    return False, "Molecule is heavy and has a protonated amine together with a carboxylate group (suggesting a zwitterion) rather than a free cation."
                elif has_positive:
                    return True, f"Heavy molecule (MW={mol_wt:.1f} Da) has net positive charge {net_charge} and a recognized positive site (though pH-dependent)."
                else:
                    return False, f"Heavy molecule has net positive charge {net_charge} but no identifiable positive site."
            else:
                return True, f"Heavy molecule (MW={mol_wt:.1f} Da) has net positive charge {net_charge} and contains a permanent cationic group."
    
    # ---------------------------
    # Case 2. net_charge < 0
    if net_charge < 0:
        return False, f"Molecule has net negative charge ({net_charge}), not a cation."
    
    # ---------------------------
    # Case 3. net_charge == 0 (potential zwitterions)
    # Count positive and negative centers
    pos_sum = sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0)
    neg_sum = -sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
    
    # If positive centers outweigh negatives, treat as cationic
    if pos_sum > neg_sum:
        return True, f"Molecule has net zero overall charge but positive centers ({pos_sum}) outweigh negative centers ({neg_sum})."
    
    # If balanced, see if a robust (permanent) cationic group is present.
    if pos_sum == neg_sum:
        if mol.HasSubstructMatch(quaternary) or mol.HasSubstructMatch(aromatic_nitrogen) or (guanidinium is not None and mol.HasSubstructMatch(guanidinium)):
            return True, "Molecule has net zero charge but contains a permanent cationic functional group."
    
    # If no rule applies, then do not classify as a cation.
    return False, f"Molecule has a net charge of {net_charge} with balanced positive ({pos_sum}) and negative ({neg_sum}) centers and no unambiguous permanent cationic substructure."

# Example usage (for testing, uncomment the following lines):
# test_smiles = [
#     # True positives:
#     "P(OCC[N+](C)(C)C)(OCC(OC(=O)CCCCCCCCCCCCC=1OC(CCCCC)=CC1C)COC(=O)CCCCCCCCCCC=2OC(=C(C2C)C)CCCCC)(O)=O",  # phosphinic acid derivative with cationic site
#     "[NH3+]CCC(=O)NCCC1=CNC=N1",  # carcininium
#     "COc1cccc(c1)[C@@]1(O)CCCC[C@@H]1C[NH+](C)C",  # (R,R)-tramadol(1+)
#     # False positives examples (should return False):
#     "OCC[NH+](C)C",  # N-dimethylethanolamine (small and pH-dependent; likely zwitterionic salt)
#     "CCCC\C=C\C\C=C\CCCCCCCC(=O)OC(CC([O-])=O)C[N+](C)(C)C"  # example for a carnitine derivative (zwitterion)
# ]
#
# for smi in test_smiles:
#     result, reason = is_cation(smi)
#     print(f"SMILES: {smi}\nResult: {result} | Reason: {reason}\n")