"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: 2-monoglyceride
A 2-monoglyceride is defined as a monoglyceride where the acyl substituent is located at position 2.
That is, the glycerol backbone (normally HOCH2–CHOH–CH2OH) has its central –OH replaced by an ester function 
so that the key substructure is: HOCH2–CH(OC(=O)R)–CH2OH.
This implementation uses a SMARTS pattern that explicitly captures the connectivity of the glycerol backbone 
with an ester at the 2‑position.
"""
from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    
    A 2-monoglyceride contains a glycerol backbone where the acyl (ester) group is attached at the 
    secondary (2-) position, i.e. HOCH2–CH(OC(=O)R)–CH2OH.
    
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
    
    # Add explicit hydrogens so that –OH groups are represented explicitly.
    mol_with_H = Chem.AddHs(mol)
    
    # Define the SMARTS pattern for a glycerol backbone that is substituted at the 2-position.
    # This pattern looks for:
    #   - a primary carbon (CH2) that bears an –OH group ([CH2;H2]([OX2H])),
    #   - connected to a secondary carbon (CH with one hydrogen, [C;H1]) that carries an ester substituent (O[C;X3](=O)[*]),
    #   - and joined to another primary CH2 group with an –OH ([CH2;H2]([OX2H])).
    pattern_smarts = "[CH2;H2]([OX2H])-[C;H1](O[C;X3](=O)[*])-[CH2;H2]([OX2H])"
    pattern = Chem.MolFromSmarts(pattern_smarts)
    if pattern is None:
        return False, "Failed to generate SMARTS pattern"
    
    # Find substructure matches for the pattern.
    matches = mol_with_H.GetSubstructMatches(pattern)
    
    if not matches:
        return False, "No 2-monoglyceride (sn-2 ester on glycerol) fragment found"
    
    # Expect exactly one glycerol backbone fragment matching the 2-monoglyceride definition.
    if len(matches) > 1:
        return False, f"Found {len(matches)} glycerol ester fragments; expected exactly one for a monoglyceride"
    
    return True, "Contains a glycerol backbone with an acyl ester at the 2-position (2-monoglyceride)"


# Example usage (for testing):
if __name__ == "__main__":
    test_smiles = [
        "O(C(CO)CO)C(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",  # MG(0:0/22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0)
        "C(\\C[C@H]1[C@H](CC([C@@H]1/C=C/[C@H](CCCCC)O)=O)O)=C\\CCCC(=O)OC(CO)CO",  # prostaglandin D2 2-glyceryl ester
        "CCCCCCCCCCCCCCCC(=O)OC(CO)CO",  # 2-palmitoylglycerol
        "CCO",  # Clearly not a 2-monoglyceride
    ]
    
    for s in test_smiles:
        valid, reason = is_2_monoglyceride(s)
        print(f"SMILES: {s}\nResult: {valid}, Reason: {reason}\n")