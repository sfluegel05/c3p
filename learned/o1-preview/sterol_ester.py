"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid nucleus pattern (cyclopentanoperhydrophenanthrene)
    # Steroid nucleus: three fused six-membered rings and one five-membered ring
    steroid_nucleus_smarts = """
    [#6]1-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-1-[#6]2-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-2-[#6]3-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-3-[#6]4-[#6]-[#6]-[#6]-[#6]-1-4
    """
    steroid_nucleus = Chem.MolFromSmarts(steroid_nucleus_smarts)
    if not mol.HasSubstructMatch(steroid_nucleus):
        return False, "No steroid nucleus found"

    # Define esterified 3-hydroxy group pattern on steroid nucleus
    # This pattern looks for an esterified oxygen at position 3 of the steroid nucleus
    esterified_3_hydroxy_smarts = """
    [#6;R1]1([#6;R1])[#6;R1][#6;R1]([#6;R1])[#6;R1][#6;R1]1
    [#6;R1]-[#8;R1]-[#6](=O)-[#6]
    """
    esterified_3_hydroxy = Chem.MolFromSmarts("""
    [#6;r6]-1-[#6;r6]-[#6;r6]-[#6;r6]-[#6;r6]-[#6;r6]-1-[#8]-[#6](=O)-[*]
    """)
    if not mol.HasSubstructMatch(esterified_3_hydroxy):
        return False, "No esterification at 3-hydroxy group of sterol nucleus found"

    return True, "Contains sterol nucleus with esterified 3-hydroxy group"