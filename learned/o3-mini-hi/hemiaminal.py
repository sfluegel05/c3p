"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: Hemi-aminal compounds.
Definition: A hemiaminal is any organic amino compound that has an amino group
and a hydroxy group attached to the same tetrahedral (sp3, non‐aromatic) carbon that 
carries exactly one hydrogen. In a correct hemiaminal motif the candidate carbon must 
have exactly three heavy-atom (non-hydrogen) neighbors: one is –OH (with the oxygen bearing 
at least one hydrogen) and one is –NH (which may be primary, secondary, or tertiary).
"""

from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule contains a hemiaminal motif.
    
    A hemiaminal here is defined as a tetrahedral (sp3, non-aromatic) carbon carrying exactly one hydrogen,
    with exactly three heavy atom (non-hydrogen) neighbors. Two of these neighbors must be:
      - An oxygen that is itself bonded to at least one hydrogen (i.e. an -OH group),
      - A nitrogen (i.e. an amino group, which can be -NH, -NHR, or -NR2).
    
    Args:
        smiles (str): SMILES representation of the molecule.
    
    Returns:
        bool: True if the molecule contains a hemiaminal motif, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to catch implicit ones.
    mol = Chem.AddHs(mol)
    
    # Iterate over all atoms and look for candidate carbons.
    for atom in mol.GetAtoms():
        # We are interested only in carbons.
        if atom.GetAtomicNum() != 6:
            continue
        
        # Must be sp3 hybridized and non-aromatic.
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or atom.GetIsAromatic():
            continue
        
        # Get the neighbors of the candidate carbon.
        neighbors = atom.GetNeighbors()
        # Count how many heavy atoms (non-H) are attached.
        heavy_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() != 1]
        # In a tetrahedral carbon with exactly one H, there must be exactly 3 heavy atom bonds.
        if len(heavy_neighbors) != 3:
            continue
        
        # Count hydrogen atoms attached to the carbon.
        # Since we added explicit hydrogens, we can simply count by atomic number.
        h_count = sum(1 for nbr in neighbors if nbr.GetAtomicNum() == 1)
        if h_count != 1:
            continue
        
        # Flags for the two required substituents.
        found_OH = False
        found_NH = False
        
        # Check each neighbor (only consider heavy atoms for the substituents).
        for nbr in heavy_neighbors:
            # Get the bond between atom and neighbor.
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            # We require that the connecting bond is a single bond.
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue
            
            # If neighbor is oxygen, check that it is part of an -OH (has at least one hydrogen bound)
            if nbr.GetAtomicNum() == 8:
                # Count hydrogens on the oxygen (explicit only now)
                h_on_o = sum(1 for sub in nbr.GetNeighbors() if sub.GetAtomicNum() == 1)
                if h_on_o >= 1:
                    found_OH = True
            
            # If neighbor is nitrogen, we consider it as the needed amino group.
            elif nbr.GetAtomicNum() == 7:
                found_NH = True
        
        if found_OH and found_NH:
            reason = ("Molecule contains a hemiaminal motif: a tetrahedral, non-aromatic sp3 carbon carrying exactly one hydrogen "
                      "and exactly three heavy-atom neighbors, with at least one substituent being a hydroxyl (-OH) group (oxygen with "
                      "≥1 hydrogen) and one being an amino (-NH, -NHR, or -NR2) group."
                      )
            return True, reason

    return False, "No hemiaminal motif (an sp3 carbon with exactly one hydrogen and exactly 3 heavy-atom neighbors, including one -OH and one -NH) was found."

# Example usage (for testing purposes):
# Uncomment the lines below to run a few examples.
# test_smiles = [
#     "NC(O)C(O)=O",          # alpha-hydroxyglycine (should be True)
#     "OC(N)CC",              # 2-Aminopropanol (should be True)
#     "C1CCCCC1"              # Cyclohexane (should be False)
# ]
# for smi in test_smiles:
#     result, explanation = is_hemiaminal(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nExplanation: {explanation}\n")