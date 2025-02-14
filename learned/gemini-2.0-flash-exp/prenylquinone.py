"""
Classifies: CHEBI:26255 prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a quinone core (5 or 6 membered ring with two carbonyls and conjugated double bonds)
    quinone_patterns = [
        Chem.MolFromSmarts("[C;$(C=O)]1[C,c](=[C,c][C,c](=[C,c]1)[C;$(C=O)])"), # 5 or 6 membered ring with two carbonyls, with substitution
        Chem.MolFromSmarts("[C;$(C=O)]1[C,c]=[C,c][C,c](=[C,c][C,c]1[C;$(C=O)])"),  # 5 or 6 membered ring with two carbonyls, with substitution
        Chem.MolFromSmarts("C1(=O)C=CC2=CC=CC=C2C1=O"), #Naphtoquinone (no substitution)
         Chem.MolFromSmarts("C1(=O)C=C[C,c]2=[C,c][C,c]=[C,c][C,c]=[C,c]2C1=O") #Naphtoquinone allowing for substition

    ]
    
    found_quinone = False
    for pattern in quinone_patterns:
      if mol.HasSubstructMatch(pattern):
        found_quinone = True
        break

    if not found_quinone:
        return False, "No quinone core found"


    # Check for at least two isoprenoid units (C-C=C-C) attached to the quinone
    # The chain can have variable lengths.
    prenyl_pattern = Chem.MolFromSmarts("[C;$(C=O)]~[C,c]!@C[CH2]C(=C)C~[CH2]C=C[CX4]")
    prenyl_matches = mol.GetSubstructMatches(prenyl_pattern)
    if len(prenyl_matches) == 0:
        return False, "No polyprenyl side chain found attached to the quinone"
    
    #Check for minimum number of carbons in the molecule to filter out false positives
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 13: #Arbitrary number based on the examples. A molecule should have at least 13 carbons to be a prenylquinone
        return False, "Too few carbons for prenylquinone"

    # Check molecular weight - prenylquinones typically > 200 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for prenylquinone"


    return True, "Contains a quinone core and a polyprenyl side chain"