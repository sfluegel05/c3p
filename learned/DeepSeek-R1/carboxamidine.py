"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    Carboxamidines have the structure RC(=NR)NR2 (a carbon double-bonded to an NR group and single-bonded to an NR2 group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the carboxamidine pattern: C connected to two nitrogens, one via double bond
    # SMARTS pattern: [CX3](=[NX2])[NX3]
    pattern = Chem.MolFromSmarts("[CX3](=[NX2])[NX3]")
    if mol.HasSubstructMatch(pattern):
        return True, "Contains carboxamidine group (RC(=NR)NR2)"
    else:
        return False, "No carboxamidine group detected"