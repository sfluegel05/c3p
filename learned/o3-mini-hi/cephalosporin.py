"""
Classifies: CHEBI:23066 cephalosporin
"""
"""
Classifies: cephalosporin
Definition: A beta‐lactam antibiotic whose core comprises a fused bicyclic ring system in which a 
4‐membered beta‐lactam ring (with one ring‐bound nitrogen and three carbons – one of which bears a carbonyl group) 
is fused to a 6‐membered dihydrothiazine ring containing at least one sulfur (or oxygen) and being non‐aromatic. 
This program uses a SMARTS query to detect a cephalosporin core.
"""
from rdkit import Chem

def is_cephalosporin(smiles: str):
    """
    Determines if a molecule is a cephalosporin based on its SMILES string.
    This function searches for a cephalosporin core defined by the SMARTS pattern:
    
        O=C1NC2SCC2C1
        
    which approximates the fused 4-membered beta-lactam and 6-membered dihydrothiazine ring system.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a cephalosporin core is detected, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS for the cephalosporin core.
    # This pattern represents a beta-lactam (4-membered) ring fused to a dihydrothiazine (6-membered) ring.
    # The pattern "O=C1NC2SCC2C1" does not include chirality or substituent details.
    cephalosporin_core_smarts = "O=C1NC2SCC2C1"
    pattern = Chem.MolFromSmarts(cephalosporin_core_smarts)
    if pattern is None:
        return False, "Error in SMARTS pattern definition"
    
    # Search for the cephalosporin core in the molecule
    if mol.HasSubstructMatch(pattern):
        return True, "Cephalosporin core structure detected (fused beta-lactam and dihydrothiazine rings found)"
    else:
        return False, "No cephalosporin core structure detected"

# Below is an example test (for demonstration purposes); in practice, one would test with the given SMILES.
if __name__ == '__main__':
    # 7beta-aminodeacetoxycephalosporanic acid (one of the known cephalosporin structures)
    test_smiles = "[H][C@]12SCC(C)=C(N1C(=O)[C@H]2N)C(O)=O"
    is_ceph, reason = is_cephalosporin(test_smiles)
    print("Is cephalosporin:", is_ceph)
    print("Reason:", reason)