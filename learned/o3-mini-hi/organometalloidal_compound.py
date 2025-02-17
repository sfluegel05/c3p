"""
Classifies: CHEBI:143084 organometalloidal compound
"""
"""
Classifies: organometalloidal compound
Definition: A compound having bonds between one or more metalloid atoms and one or more carbon atoms 
             of an organyl group.
"""
from rdkit import Chem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    The compound must have at least one bond between a carbon atom and a metalloid (e.g., B, Si, Ge, As, Sb, Te).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a list of atomic numbers for metalloids
    # Boron (5), Silicon (14), Germanium (32), Arsenic (33), Antimony (51), Tellurium (52)
    metalloids = {5, 14, 32, 33, 51, 52}
    
    # Iterate over all bonds in the molecule.
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        a1_atomic_num = atom1.GetAtomicNum()
        a2_atomic_num = atom2.GetAtomicNum()
        
        # Check if one atom is carbon (atomic number 6) and the other is a metalloid.
        if (a1_atomic_num == 6 and a2_atomic_num in metalloids) or (a2_atomic_num == 6 and a1_atomic_num in metalloids):
            # Found an organometalloidal bond.
            return True, f"Found bond between carbon (atomic number 6) and metalloid (atomic number {a2_atomic_num if a1_atomic_num == 6 else a1_atomic_num})."
    
    # If no such bond is found, the compound does not fulfill the definition.
    return False, "No bond between carbon and a metalloid atom found."

# Example usage:
if __name__ == "__main__":
    example_smiles = "[As](=O)(C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCC)(C)C"
    result, reason = is_organometalloidal_compound(example_smiles)
    print(f"Result: {result}\nReason: {reason}")