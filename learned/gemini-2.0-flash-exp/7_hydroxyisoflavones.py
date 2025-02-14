"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone has a hydroxyl group at the 7-position of the isoflavone core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise.
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for the 7-hydroxyisoflavone core structure
    # using generic atoms (a), allowing for substitutions
    # 'c1cc(O)cc2c(c1)c(=O)oc2' or 'c1cc(O)cc2c(c1)oc(=O)c2' with substitutents
    
    core_pattern1 = Chem.MolFromSmarts('a1c(O)cc(a)c2c(a1)c(=O)oc2')
    core_pattern2 = Chem.MolFromSmarts('a1c(O)cc(a)c2c(a1)oc(=O)c2')
    
    if not (mol.HasSubstructMatch(core_pattern1) or mol.HasSubstructMatch(core_pattern2)):
        return False, "Molecule does not have the 7-hydroxyisoflavone core structure"

    # Count the number of hydroxy groups on the benzene ring of the core
    benzene_ring_pattern = Chem.MolFromSmarts('c1cc(O)cc2')
    benzene_matches = mol.GetSubstructMatches(benzene_ring_pattern)
    
    if benzene_matches:
        for match in benzene_matches:
           
            hydroxy_count = 0
            for atom_index in match:
                atom = mol.GetAtomWithIdx(atom_index)
                if atom.GetSymbol() == 'O' and atom.GetTotalValence() == 1 :
                        hydroxy_count += 1

            if hydroxy_count != 1:
                return False, "The benzene ring must have only 1 hydroxy group in position 7"
    else:
        return False, "Isoflavone core must have the benzene with at least one hydroxyl group"


    return True, "Molecule has the 7-hydroxyisoflavone core with a hydroxyl group at the 7-position"