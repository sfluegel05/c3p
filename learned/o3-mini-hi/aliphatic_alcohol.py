"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
"""
Classifies: An aliphatic alcohol, defined as 'An alcohol derived from an aliphatic compound.'
This function checks if the molecule (given as a SMILES string) contains at least one -OH group 
attached to an sp3, non-aromatic carbon, which defines an aliphatic alcohol group.
"""

from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is defined here as one containing at least one hydroxyl (-OH) group 
    attached to a saturated, non-aromatic (sp3) carbon.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as an aliphatic alcohol, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms; look for oxygen atoms with at least one hydrogen attached (i.e., -OH group)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            # Check if the atom has at least one hydrogen (either explicit or implicit)
            # GetTotalNumHs() returns the total number of attached hydrogens.
            if atom.GetTotalNumHs() > 0:
                # Check neighbor atoms that might be the carbon to which the hydroxyl is attached.
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        # Check if the carbon is sp3 (saturated) and not aromatic
                        if neighbor.GetHybridization() == Chem.rdchem.HybridizationType.SP3 and not neighbor.GetIsAromatic():
                            return True, "Found -OH group attached to an aliphatic (sp3, non-aromatic) carbon"
    
    return False, "No aliphatic -OH group found; either there is no hydroxyl group or it is not attached to a saturated carbon."
    
# Below are some tests (not required) to verify our function with example SMILES.
if __name__ == "__main__":
    test_smiles = [
        "O=C1OC([C@@H](O)\\C=C/C=C/C)CC1",  # Sapinofuranone A
        "CCCCCCC(C)O",                    # octan-2-ol
        "CCCCCCCCCCCCCCCCCCCCCCO",         # tricosan-1-ol
        "OC1=CC=CC=C1"                    # phenol (should fail: aromatic alcohol)
    ]
    
    for sm in test_smiles:
        valid, msg = is_aliphatic_alcohol(sm)
        print(f"SMILES: {sm}\nClassification: {valid}\nReason: {msg}\n")