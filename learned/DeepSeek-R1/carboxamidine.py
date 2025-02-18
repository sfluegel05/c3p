"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

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

    # SMARTS pattern for carboxamidine group
    # [CX3;!a]: non-aromatic carbon with three bonds
    # =[NX2;D2]: double bond to N with two connections (C and R/H)
    # [NX3;D3]: single bond to N with three connections (C and two R/H)
    pattern = Chem.MolFromSmarts("[CX3;!a](=[NX2;D2])[NX3;D3]")
    
    if mol.HasSubstructMatch(pattern):
        return True, "Contains carboxamidine group (RC(=NR)NR2)"
    else:
        return False, "No carboxamidine group detected"