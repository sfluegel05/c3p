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
      - an oxygen that is bound as -OH (i.e. oxygen that carries at least one hydrogen), and 
      - a nitrogen (i.e. an amino group, which can be primary, secondary, or tertiary).

    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule contains a hemiaminal motif, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES and add explicit hydrogens for proper counting.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Iterate over atoms looking for candidate carbon atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # only consider carbon atoms
        
        # Candidate carbon must be sp3 hybridized and non-aromatic.
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3 or atom.GetIsAromatic():
            continue
        
        # Using explicit hydrogens, the total degree of the carbon should be 4.
        if atom.GetDegree() != 4:
            continue
        
        # Count hydrogens attached to carbon using GetTotalNumHs.
        h_count = atom.GetTotalNumHs()
        if h_count != 1:
            continue  # must have exactly one hydrogen
        
        # Count heavy (non-hydrogen) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 3:
            continue

        # Ensure all bonds from candidate carbon are single bonds.
        bonds_all_single = all(mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()).GetBondType() 
                               == Chem.rdchem.BondType.SINGLE for nbr in heavy_neighbors)
        if not bonds_all_single:
            continue
            
        # Look for the required substituents: one oxygen that is -OH and one nitrogen.
        found_OH = False
        found_N = False
        for nbr in heavy_neighbors:
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                continue  # only consider single bonds
            
            # If neighbor is oxygen, check that it is part of an -OH (at least one hydrogen attached).
            if nbr.GetAtomicNum() == 8:
                # Use GetTotalNumHs on oxygen to determine attached hydrogen count.
                if nbr.GetTotalNumHs() >= 1:
                    found_OH = True
            # If neighbor is nitrogen, flag this substituent.
            elif nbr.GetAtomicNum() == 7:
                found_N = True
                
        if found_OH and found_N:
            reason = ("Molecule contains a hemiaminal motif: a tetrahedral, non-aromatic sp3 carbon carrying exactly one hydrogen "
                      "and exactly three heavy-atom neighbors, with one substituent being a hydroxyl (-OH) group (oxygen bearing at least one hydrogen) "
                      "and one being an amino (-NH, -NHR, or -NR2) group.")
            return True, reason

    return False, ("No hemiaminal motif (an sp3 carbon with exactly one hydrogen and exactly 3 heavy-atom neighbors, "
                   "including one -OH substituent and one -NH substituent) was found.")

# Example usage (uncomment to test):
# test_smiles = [
#     "NC(O)C(O)=O",           # alpha-hydroxyglycine (should be True)
#     "OC(N)CC",               # 2-Aminopropanol (should be True)
#     "C1CCCCC1"               # Cyclohexane (should be False)
# ]
# for smi in test_smiles:
#     result, explanation = is_hemiaminal(smi)
#     print(f"SMILES: {smi}\nResult: {result}\nExplanation: {explanation}\n")