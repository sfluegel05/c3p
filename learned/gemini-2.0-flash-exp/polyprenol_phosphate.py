"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is a polyprenol chain esterified with a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str). True if molecule is a polyprenol phosphate, False otherwise, along with a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for phosphate group (allowing for deprotonation)
    phosphate_pattern = Chem.MolFromSmarts("[P](=[O])([O])([O,H])[O,H]") #P with 3 or 4 oxygens.
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Define SMARTS for isoprenoid unit within polyprenol chain
    isoprenoid_unit = Chem.MolFromSmarts("C[C](C)=C[C]") #one methyl branch
    isoprenoid_matches = mol.GetSubstructMatches(isoprenoid_unit)
    
    if len(isoprenoid_matches) < 1:
         return False, "Not enough isoprenoid units in chain"

    #Check for connection between phosphate and the chain
    connectivity_pattern = Chem.MolFromSmarts("[CX4,CX3]O[P]")
    connectivity_matches = mol.GetSubstructMatches(connectivity_pattern)
    
    if len(connectivity_matches) < 1 :
         return False, "Phosphate group is not connected to chain via an oxygen"
    
    
    # Check chain length - count carbons and make sure it's a reasonable multiple of 5 (length of isoprene)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, "Too few carbons for polyprenol"

    return True, "Contains a polyprenol chain esterified with a phosphate group"