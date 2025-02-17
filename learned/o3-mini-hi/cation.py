"""
Classifies: CHEBI:36916 cation
"""
"""
Classifies chemical entities of the class cation:
A monoatomic or polyatomic species having one or more elementary charges of the proton.
This improved classifier uses overall net formal charge along with SMARTS matching
for permanent (non-pH‐dependent) cationic functional groups.
"""

from rdkit import Chem

def is_cation(smiles: str):
    """
    Determines if a molecule is a cation based on its SMILES string.
    
    A molecule is classified as a cation if:
      - It has a net formal charge greater than zero and at least one cationic site is recognized,
    OR
      - It shows net zero overall formal charge but contains a permanent cationic functional group
        (such as quaternary ammonium or aromatic nitrogen) or the sum of positive formal charges
        outweighs negative ones.
        
    Note:
      Because many molecules are zwitterions, the algorithm first checks the net charge.
      Then for species with net zero, it compares the total positive versus negative centers.
      It also uses refined SMARTS patterns: for example, protonated primary amines (e.g. [NH3+])
      may be pH dependent; so only if such groups are present without a typical carboxylate (C(=O)[O-])
      will they be taken as evidence of a standalone cation.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a cation, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute overall net formal charge
    net_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    
    # Define several SMARTS patterns for (potentially) permanent cationic groups
    # Here the quaternary ammonium (no attached hydrogen) is generally pH‐independent:
    quaternary_ammonium = Chem.MolFromSmarts("[N+](C)(C)(C)")
    # Aromatic or heterocyclic nitrogen positive centers (e.g. thiazolium, pyridinium)
    aromatic_nitrogen = Chem.MolFromSmarts("[n+]")
    # Protonated amines (often pH-dependent); these appear in many zwitterions.
    protonated_amine = Chem.MolFromSmarts("[NH3+]")
    
    # Also define anionic groups that are common counter ions:
    carboxylate = Chem.MolFromSmarts("C(=O)[O-]")
    phosphate = Chem.MolFromSmarts("P(=O)([O-])")
    
    # ---------------------------
    # Case 1. net_charge > 0:
    if net_charge > 0:
        # When overall net charge is positive, require that a recognizable cationic group is present.
        if (mol.HasSubstructMatch(quaternary_ammonium) or 
            mol.HasSubstructMatch(aromatic_nitrogen) or 
            mol.HasSubstructMatch(protonated_amine)):
            return True, f"Molecule has a net positive charge of {net_charge} and contains a cationic site."
        else:
            return False, f"Molecule has a net positive charge of {net_charge} but lacks a recognizable permanent cationic group."
    
    # ---------------------------
    # Case 2. net_charge < 0:
    if net_charge < 0:
        return False, f"Molecule has a net negative charge of {net_charge} (not a cation)."
    
    # ---------------------------
    # Case 3. net_charge == 0 (zwitterionic candidates)
    # Sum of positive and negative formal charges (as absolute amounts)
    pos_sum = sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0)
    neg_sum = -sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
    
    # If the positive centers outnumber the negative ones, we consider it cationic.
    if pos_sum > neg_sum:
        return True, f"Molecule has a net zero charge but positive centers ({pos_sum}) outweigh negative centers ({neg_sum})."
    
    # If the totals balance, look for a definitive permanent cationic substructure.
    if pos_sum == neg_sum:
        # Quaternary ammonium or aromatic nitrogen are less pH‐dependent indicators.
        if mol.HasSubstructMatch(quaternary_ammonium) or mol.HasSubstructMatch(aromatic_nitrogen):
            return True, "Molecule contains a permanent cationic functional group (e.g., quaternary ammonium or aromatic nitrogen) despite net zero charge."
        # Next, if a protonated amine is present and no carboxylate is found (which could indicate an acid–base zwitterion), classify as cation.
        if mol.HasSubstructMatch(protonated_amine) and not mol.HasSubstructMatch(carboxylate):
            return True, "Molecule contains a protonated amine without a carboxylate group, indicative of a cation despite net zero charge."
    
    # If none of the above conditions are met, we do not classify the entity as a cation.
    return False, f"Molecule has a net charge of {net_charge} and no identified permanent cationic group or unbalanced charge centers."

# Example usage (uncomment to test):
# test_smiles = [
#     # Examples that should be classified as cations:
#     "P(OCC[N+](C)(C)C)(OCC(OC(=O)CCCCCCCCCCCCC=1OC(CCCCC)=CC1C)COC(=O)CCCCCCCCCCC=2OC(=C(C2C)C)CCCCC)(O)=O",  # phosphinic acid derivative (cationic functional group)
#     "[NH3+]CCC(=O)NCCC1=CNC=N1",  # carcininium
#     "COc1cccc(c1)[C@@]1(O)CCCC[C@@H]1C[NH+](C)C",  # (R,R)-tramadol(1+)
#     # Examples that should not be classified as cations:
#     "O=C(NCC[NH3+])C[C@](CC(=O)NC[C@@H](C([O-])=O)[NH3+])(C([O-])=O)O",  # zwitterionic peptide-like species
#     "OCC[NH+](C)C"  # N-dimethylethanolamine (a false positive case)
# ]
#
# for smi in test_smiles:
#     result, reason = is_cation(smi)
#     print(f"SMILES: {smi}\nResult: {result} | Reason: {reason}\n")