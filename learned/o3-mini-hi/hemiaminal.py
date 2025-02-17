"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: Hemi-aminal compounds
Definition: Any organic amino compound that has an amino group and a hydroxy group 
attached to the same tetrahedral (sp3) carbon that still carries one hydrogen.
Hemiaminals are typically intermediates in imine formation.
This code uses connectivity and hydrogen‐count checks to improve discrimination.
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule contains a hemiaminal motif.
    A hemiaminal is defined here as a tetrahedral (sp3) carbon, carrying exactly one 
    hydrogen (so that it is not an aminal), that is connected via single bonds to both 
    a hydroxyl group (-OH) and an amino group (-NH, -NHR, or -NR2).
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule contains a hemiaminal motif, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    hemiaminal_found = False
    reason_info = ""
    
    # Loop over all carbon atoms in the molecule. We restrict to sp3 carbons.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # For a typical hemiaminal carbon, we expect it to be tetrahedral with exactly three heavy-atom neighbors
        # (so that adding implicit hydrogens gives a total of 4 bonds) and exactly one hydrogen.
        if (atom.GetTotalNumHs() != 1) or (len(atom.GetNeighbors()) != 3):
            continue
        
        # Initialize flags for finding the required substituents.
        has_hydroxyl = False
        has_amine = False
        
        # Examine bonds from this carbon.
        for neighbor in atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            # We only want single bonds.
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            
            # Check for a hydroxyl group: oxygen that carries at least one hydrogen.
            if neighbor.GetAtomicNum() == 8:
                # In RDKit hydrogens are often implicit; get the total H count.
                if neighbor.GetTotalNumHs() >= 1:
                    has_hydroxyl = True
            # Check for an amine group: nitrogen.
            elif neighbor.GetAtomicNum() == 7:
                has_amine = True
        
        # If both -OH and -N substituents are found, we consider this a hemiaminal carbon.
        if has_hydroxyl and has_amine:
            hemiaminal_found = True
            reason_info = ("Molecule contains a hemiaminal motif: a tetrahedral sp3 carbon (with one hydrogen) "
                           "directly bound to both a hydroxyl group (-OH) and an amino group (-NH, -NHR, or -NR2).")
            break  # stop at first found match
    
    if hemiaminal_found:
        return True, reason_info
    else:
        return False, "No hemiaminal motif (a tetrahedral carbon bearing both -OH and -N with one hydrogen) was found."

# Example usage (for testing, uncomment to run):
# test_smiles = "NC(O)C(O)=O"  # Example: alpha-hydroxyglycine
# result, explanation = is_hemiaminal(test_smiles)
# print(result, explanation)