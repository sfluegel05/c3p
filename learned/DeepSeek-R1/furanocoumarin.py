"""
Classifies: CHEBI:24128 furanocoumarin
"""
"""
Classifies: CHEBI:3616 furanocoumarin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin consists of a furan ring fused with a coumarin (benzene fused to a pyrone ring).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define coumarin substructure (benzene fused to pyrone)
    coumarin_pattern = Chem.MolFromSmarts('[#6]1[#6][#6][#6][#6][#6]1-[#6]2[#6](=[#8])-[#8][#6]=[#6]2')
    if not mol.HasSubstructMatch(coumarin_pattern):
        return False, "No coumarin substructure found"

    # Define furan substructure (o1cccc1)
    furan_pattern = Chem.MolFromSmarts('[#8]1-,=[#6]-[#6]-,=[#6]-[#6]-,=1')
    furan_matches = mol.GetSubstructMatches(furan_pattern)
    if not furan_matches:
        return False, "No furan ring found"

    # Check if any furan is fused to the coumarin
    coumarin_matches = mol.GetSubstructMatches(coumarin_pattern)
    for cm in coumarin_matches:
        cm_atoms = set(cm)
        for fm in furan_matches:
            # Check adjacent atoms in furan that are part of coumarin
            furan_submol = Chem.MolFromSmarts('[#8]1-,=[#6]-[#6]-,=[#6]-[#6]-,=1')
            furan_bonds = furan_submol.GetBonds()
            for bond in furan_bonds:
                begin_idx = bond.GetBeginAtomIdx()
                end_idx = bond.GetEndAtomIdx()
                actual_begin = fm[begin_idx]
                actual_end = fm[end_idx]
                if actual_begin in cm_atoms and actual_end in cm_atoms:
                    return True, "Furan ring fused to coumarin structure"

    return False, "Furan not fused to coumarin"