"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
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

    # Define SMARTS pattern for phosphate group (allowing for mono and di-phosphate)
    phosphate_pattern = Chem.MolFromSmarts("[P](=[O])([O])([O,H])[O,H]")
    diphosphate_pattern = Chem.MolFromSmarts("[P](=[O])([O])([O,H])[O][P](=[O])([O])([O,H])[O,H]")
    if not mol.HasSubstructMatch(phosphate_pattern) and not mol.HasSubstructMatch(diphosphate_pattern):
       return False, "No phosphate group found"

    # Define SMARTS for a repeating isoprene unit, with a minimum chain length
    isoprene_pattern = Chem.MolFromSmarts("C(C)=[C]C[C]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2: #minimum of two isoprene units
        return False, "Not enough repeating isoprene units in chain"

    #Check connectivity, look for -O-P. At least one connection.
    connectivity_pattern = Chem.MolFromSmarts("[OX2][P]")
    connectivity_matches = mol.GetSubstructMatches(connectivity_pattern)
    if len(connectivity_matches) < 1:
       return False, "Phosphate group is not connected to chain via oxygen"
    

    # Check number of oxygen atoms - should not be a carbohydrate
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 10:
        # Check for simple carbohydrate, e.g., [C][C]([O])[C]([O])[C]([O])[C]([O])[C][O]
        carb_pattern = Chem.MolFromSmarts("[CX4][CX4]([OX2])[CX4]([OX2])[CX4]([OX2])[CX4]([OX2])[CX4][OX2]")
        if mol.HasSubstructMatch(carb_pattern):
           return False, "Too many oxygen atoms, likely a carbohydrate"

    return True, "Contains a polyprenol chain esterified with a phosphate group"