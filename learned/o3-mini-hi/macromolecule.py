"""
Classifies: CHEBI:33839 macromolecule
"""
#!/usr/bin/env python3
"""
Classifies: Macromolecule

Definition:
  "A macromolecule is a molecule of high relative molecular mass, 
   the structure of which essentially comprises the multiple repetition 
   of units derived, actually or conceptually, from molecules of low relative molecular mass."

This improved heuristic makes use of multiple thresholds:
  • If the molecular weight is very high (>=1000 Da), we classify the molecule as macromolecular.
  • If weight is between 800 and 1000 Da, we also require at least 4 amide bonds and a reasonable atom count (>=40).
  • If weight is between 500 and 800 Da, we require at least 5 amide bonds and an atom count >=50.
Note: This is still an approximate method and may miss some macromolecules or mis‐classify others.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    
    Heuristics used:
      1. If molecular weight >= 1000 Da, the molecule is macromolecular.
      2. If molecular weight is between 800 and 1000 Da, then there must be at least 4 amide bonds
         and at least 40 atoms.
      3. If molecular weight is between 500 and 800 Da, then there must be at least 5 amide bonds
         and at least 50 atoms.
         
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if molecule is classified as a macromolecule, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    
    # Count total number of atoms (a rough proxy for size).
    num_atoms = mol.GetNumAtoms()
    
    # Count number of amide bonds using a SMARTS search.
    # We use a pattern for a C(=O)N group.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amide = len(amide_matches)
    
    # Apply our revised thresholds.
    if mol_wt >= 1000:
        reason = f"Molecular weight is {mol_wt:.1f} Da which is high for a macromolecule."
        if num_amide > 0:
            reason += f" Found {num_amide} amide bonds."
        return True, reason
    elif 800 <= mol_wt < 1000:
        if num_amide >= 4 and num_atoms >= 40:
            reason = (f"Molecular weight is {mol_wt:.1f} Da, and the molecule has {num_atoms} atoms "
                      f"with {num_amide} amide bonds, suggesting repeated units.")
            return True, reason
        else:
            reason = (f"Molecular weight is {mol_wt:.1f} Da but only {num_amide} amide bonds "
                      f"and {num_atoms} atoms were found; stricter criteria for macro repetition are not met.")
            return False, reason
    elif 500 <= mol_wt < 800:
        if num_amide >= 5 and num_atoms >= 50:
            reason = (f"Molecular weight is {mol_wt:.1f} Da, and the molecule has {num_atoms} atoms "
                      f"with {num_amide} amide bonds, suggesting multiple repeat units.")
            return True, reason
        else:
            reason = (f"Molecular weight is {mol_wt:.1f} Da with {num_atoms} atoms and only {num_amide} amide bonds; "
                      "criteria for a macromolecule are not met.")
            return False, reason
    else:
        reason = (f"Molecular weight is {mol_wt:.1f} Da with {num_atoms} atoms and only {num_amide} amide bonds; "
                  "criteria for macromolecule (high MW and repeated units) are not met.")
        return False, reason

# Example usage:
if __name__ == '__main__':
    # Test with example SMILES. 
    # (This is a short peptide that, under the new criteria, is unlikely to be macromolecular.)
    test_smiles = "CC(C)C(=O)NC(Cc1ccccc1)C(=O)NC(CO)C(=O)NCC(=O)O"
    flag, explanation = is_macromolecule(test_smiles)
    print(flag)
    print(explanation)