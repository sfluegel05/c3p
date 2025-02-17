"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: Amine
Definition: A compound formally derived from ammonia by replacing one, two or three hydrogen atoms
by hydrocarbyl groups.
Examples of valid amines include primary, secondary, and tertiary amines,
provided that the nitrogen is not directly bound to a carbonyl carbon (i.e. not an amide).
"""

from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule contains an amine functional group based on its SMILES string.
    
    A valid amine here is defined as a nitrogen atom (atomic number 7) that is derived from ammonia
    (i.e. has three substituents overall, whether by bonds or implicit hydrogens) and is not directly
    bonded to a carbonyl carbon (which would indicate an amide functionality).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if at least one qualifying amine group is found, False otherwise.
        str: Explanation or reason for the classification.
    """
    
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # skip non-nitrogen atoms
        
        # Count substituents: sum of explicit neighbors (degree) and implicit hydrogens.
        total_substituents = atom.GetDegree() + atom.GetTotalNumHs()
        # We require that the nitrogen (if derived from ammonia) has exactly 3 substituents.
        # (Primary: 1 bond + 2 H, Secondary: 2 bonds + 1 H, Tertiary: 3 bonds + 0 H)
        if total_substituents != 3:
            continue
        
        # Exclude nitrogens bound to a carbonyl carbon (indicative of an amide functionality)
        skip_atom = False
        for neighbor in atom.GetNeighbors():
            # Check if the neighbor is a carbon that might be part of a carbonyl.
            if neighbor.GetAtomicNum() == 6:
                for bond in neighbor.GetBonds():
                    # Check if bond is a double bond (order 2) and the other atom is oxygen.
                    if bond.GetBondTypeAsDouble() == 2.0:
                        other_atom = bond.GetOtherAtom(neighbor)
                        if other_atom.GetAtomicNum() == 8:
                            skip_atom = True
                            break
                if skip_atom:
                    break
        if skip_atom:
            continue
        
        # If we have found a nitrogen which is not in an amide and has the proper number of substituents, we qualify.
        return True, "Found at least one nitrogen derived from ammonia (with three substituents) that is not an amide."
    
    # If no nitrogen qualifies as an amine
    return False, "No qualifying amine functional group found."

# The code above can be directly tested with the provided examples.
# For example:
# result, reason = is_amine("CNc1ccccc1")
# print(result, reason)