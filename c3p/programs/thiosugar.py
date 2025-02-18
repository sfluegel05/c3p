"""
Classifies: CHEBI:73754 thiosugar
"""
"""
Classifies: a thiosugar
Definition: A carbohydrate derivative in which one or more of the oxygens or hydroxy 
groups of the parent carbohydrate is replaced by sulfur or -SR, where R can be hydrogen 
or any group.
"""

from rdkit import Chem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    The function first searches for sugar ring motifs (pyranose and furanose rings) – these
    are defined by rings containing one heteroatom (usually oxygen) and the remaining carbons. 
    It then examines whether any of the positions that normally are oxygen or hydroxyl groups have been
    replaced by sulfur.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a thiosugar, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string. Return if invalid.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for typical sugar rings:
    # Pyranose: six-membered ring with one heteroatom in the ring.
    pyranose_smarts = "[O;R]1[C;R][C;R][C;R][C;R][C;R]1"
    # Furanose: five-membered ring with one heteroatom.
    furanose_smarts = "[O;R]1[C;R][C;R][C;R][C;R]1"
    
    pyranose = Chem.MolFromSmarts(pyranose_smarts)
    furanose = Chem.MolFromSmarts(furanose_smarts)
    
    # Find sugar ring matches (could be multiple saccharide moieties)
    matches = mol.GetSubstructMatches(pyranose) + mol.GetSubstructMatches(furanose)
    if not matches:
        return False, "No sugar ring structure found in molecule"
    
    # For each sugar ring candidate, check for thio substitution:
    for match in matches:
        # 'match' is a tuple of atom indices that are part of the ring.
        # Flag if we find S in the ring or as substituent on a ring carbon.
        # First, check the ring atoms themselves.
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() == 'S':
                return True, "Sugar ring contains sulfur in place of oxygen"
        
        # Next, check the non-ring substituents from ring carbon atoms.
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            # We expect typical sugars to have hydroxyl substituents on carbons.
            # If any neighbor (not in the ring) is sulfur, this indicates thio substitution.
            if atom.GetSymbol() == 'C':
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in match and nbr.GetSymbol() == 'S':
                        return True, "Sugar ring has thio substitution at a hydroxyl-bearing carbon"
    # If none of the sugar rings have a sulfur replacement.
    return False, "Sugar rings detected but no thio substitution found"
    
# Example test cases (uncomment to test)
# examples = {
#     "butylglucosinolic acid": "[C@H]1(O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO)S/C(=N\\OS(O)(=O)=O)/CCCC",
#     "celesticetin": "COC(C)C(NC(=O)[C@@H]1CCCN1C)[C@H]2O[C@H](SCCOC(=O)C3=CC=CC=C3O)[C@H](O)[C@@H](O)[C@H]2O",
#     "6-thio-beta-D-galactose": "[C@@H]1([C@@H]([C@H]([C@H]([C@H](O1)CS)O)O)O)O"
# }
# for name, smi in examples.items():
#     is_thio, reason = is_thiosugar(smi)
#     print(f"{name}: {is_thio} ({reason})")