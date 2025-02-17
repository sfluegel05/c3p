"""
Classifies: CHEBI:33839 macromolecule
"""
#!/usr/bin/env python3
"""
Classifies: Macromolecule
Definition: "A macromolecule is a molecule of high relative molecular mass, the structure
of which essentially comprises the multiple repetition of units derived, actually or conceptually,
from molecules of low relative molecular mass."
This heuristic checks the molecular weight and looks for multiple repeating amide bonds.
Note: This is an approximate method. The program may not capture all macromolecules or might mis‐classify some.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    We use two heuristics:
      1. A very high molecular weight (e.g. >= 800 Da) is typical of macromolecules.
      2. Even if the molecular weight is moderate (>=500 Da) if there are multiple (>=3) amide bonds,
         this suggests repeated units (commonly found in peptides/polymers).
    
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
    
    # Count number of atoms (a rough measure of size).
    num_atoms = mol.GetNumAtoms()
    
    # Define a SMARTS pattern for an amide bond.
    # Note: This pattern looks for a C=O attached to an N.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    num_amide = len(amide_matches)
    
    # Heuristic classification.
    # If molecular weight is very high, then we classify as macromolecule.
    if mol_wt >= 800:
        reason = f"Molecular weight is {mol_wt:.1f} Da which is high for a macromolecule."
        if num_amide >= 2:
            reason += f" Also, found {num_amide} amide bonds suggesting repeated units."
        return True, reason
    # Even if the molecular weight is a bit lower,
    # if there are several amide bonds (indicating repeating units) and a modest number of atoms, 
    # we may still classify it as a macromolecule.
    elif mol_wt >= 500 and num_amide >= 3 and num_atoms >= 20:
        reason = (f"Molecular weight is {mol_wt:.1f} Da, which is moderate, but found {num_amide} amide bonds "
                  f"in a molecule with {num_atoms} atoms, suggesting multiple repeat units.")
        return True, reason
    else:
        reason = (f"Molecular weight is {mol_wt:.1f} Da with {num_atoms} atoms and only {num_amide} amide bonds; "
                  "criteria for macromolecule (high MW and repeating units) are not met.")
        return False, reason

# Example usage:
if __name__ == '__main__':
    # Example SMILES for a peptide macromolecule (this is just an example; real macromolecule SMILES can be very long)
    test_smiles = "CC(C)C(=O)NC(Cc1ccccc1)C(=O)NC(CO)C(=O)NCC(=O)O"  # A short peptide (not really macromolecular)
    flag, explanation = is_macromolecule(test_smiles)
    print(flag)
    print(explanation)