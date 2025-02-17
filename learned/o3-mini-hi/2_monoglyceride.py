"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: 2-monoglyceride
A 2-monoglyceride is defined as a monoglyceride in which the acyl substituent is located at position 2.
That is, the glycerol backbone (normally HOCH2–CHOH–CH2OH) has its central hydroxyl replaced by an ester function so 
that the key substructure is: CH2–OH – CH(OC(=O)R) – CH2–OH.
This program uses an improved SMARTS pattern (and adds explicit hydrogens) 
to better detect the expected fragment in positive examples.
"""

from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    
    A 2-monoglyceride contains a glycerol backbone where the acyl (ester) group is attached 
    at the secondary (2-) position (i.e., HOCH2–CH(OC(=O)R)–CH2OH).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 2-monoglyceride, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that –OH groups are visible to substructure matching.
    mol_with_H = Chem.AddHs(mol)
    
    # Define a SMARTS pattern representing the glycerol backbone with an ester at the 2-position.
    # The pattern looks for a fragment with:
    #  - a primary carbon (CH2) with two hydrogens bound to an –OH group,
    #  - a secondary carbon (CH with one hydrogen) attached to an oxygen (which is in turn attached to a carbonyl),
    #  - and another primary CH2–OH.
    # The acyl (ester) part is represented as O[C;X3](=O)[*], where [*] is any substituent.
    pattern_smarts = "[CH2;H2][OH]-[C;H1](O[C;X3](=O)[*])-[CH2;H2][OH]"
    pattern = Chem.MolFromSmarts(pattern_smarts)
    if pattern is None:
        return False, "Failed to generate SMARTS pattern"
    
    # Find all substructure matches for the pattern.
    matches = mol_with_H.GetSubstructMatches(pattern)
    if not matches:
        return False, "No 2-monoglyceride (sn-2 ester on glycerol) fragment found"
    
    # For a proper monoglyceride we expect exactly one matching glycerol (sn-2 ester) fragment.
    if len(matches) > 1:
        return False, f"Found {len(matches)} glycerol ester fragments; expected exactly one for a monoglyceride"
    
    return True, "Contains a glycerol backbone with an acyl ester at the 2-position (2-monoglyceride)"


# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        "O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",  # MG(0:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0)
        "CCCCCCCCCCCCCCCC(=O)OC(CO)CO",  # 2-palmitoylglycerol
        "O(C(=O)CCCCCCCCCCCCCCCCCCCCC)C(CO)CO",  # MG(0:0/22:0/0:0)
        "CCO",  # obviously not a 2-monoglyceride
    ]
    
    for s in test_smiles:
        valid, reason = is_2_monoglyceride(s)
        print(f"SMILES: {s}\nResult: {valid}, Reason: {reason}\n")