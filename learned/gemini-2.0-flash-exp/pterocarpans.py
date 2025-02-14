"""
Classifies: CHEBI:26377 pterocarpans
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    A pterocarpan has a 6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton.
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the core 6H-[1]benzofuro[3,2-c]chromene skeleton, without stereochemistry constraints
    core_pattern = Chem.MolFromSmarts('C12Oc3c4ccccc4Oc4c3cc5ccccc5C21')

    if core_pattern is None:
       return False, "Invalid SMARTS string"

    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core 6H-[1]benzofuro[3,2-c]chromene skeleton not found."
    
    # Verify 6a and 11a positions are dihydro
    match = mol.GetSubstructMatch(core_pattern)
    
    atom_6a = mol.GetAtomWithIdx(match[0])
    atom_11a = mol.GetAtomWithIdx(match[-1])

    if (atom_6a.GetTotalNumHs() != 1 or atom_11a.GetTotalNumHs() != 1):
        return False, "6a and 11a are not dihydro"
    

    #check the presence of two fused benzene rings
    benzene_pattern = Chem.MolFromSmarts('c1ccccc1')
    benzene_matches = mol.GetSubstructMatches(benzene_pattern)

    if len(benzene_matches) < 2:
         return False, f"Found {len(benzene_matches)} benzene rings, needs at least 2"

    #Check the presence of the furan ring
    furan_pattern = Chem.MolFromSmarts('c1occc1')
    if not mol.HasSubstructMatch(furan_pattern):
       return False, "No furan ring found."
    
    return True, "Molecule contains the core 6H-[1]benzofuro[3,2-c]chromene skeleton with dihydro 6a and 11a positions"