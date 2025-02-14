"""
Classifies: CHEBI:25409 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a monoterpenoid, False otherwise
         str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 130 or mol_wt > 350: # Adjusted range
        return False, f"Molecular weight {mol_wt} is not within the typical range of monoterpenoids (130-350)."

    # Check carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8 or c_count > 15: # Adjusted range to account for some modifications
        return False, f"Too few or too many carbons: {c_count}. Monoterpenoids typically have around 10 carbons."
    
    # Check O:C ratio
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if c_count > 0 and (o_count / c_count) > 0.5:
        return False, "Monoterpenoids usually have a low ratio of oxygen to carbon atoms."

    # Common monoterpenoid substructures
    menthane_pattern = Chem.MolFromSmarts("[C]1([C])[C]([C])([C])[C]([C])[C]1")
    pinane_pattern = Chem.MolFromSmarts("C1[C]2[C]([C]1([C])C)[C]([C])([C])CC2")
    bornane_pattern = Chem.MolFromSmarts("C1[C]2[C]([C]1([C])C)[C]([C])([C])CC2")
    thujane_pattern = Chem.MolFromSmarts("C1[C]2[C]([C]1([C])C)[C]([C])([C])CC2")
    p_menthene_pattern1 = Chem.MolFromSmarts("C[C]1([C])[C]([C])([C])[C]([C])=[C]1")
    p_menthene_pattern2 = Chem.MolFromSmarts("C[C]1([C])=[C]([C])[C]([C])[C]([C])[C]1")

    if (mol.HasSubstructMatch(menthane_pattern) or
            mol.HasSubstructMatch(pinane_pattern) or
            mol.HasSubstructMatch(bornane_pattern) or
            mol.HasSubstructMatch(thujane_pattern) or
            mol.HasSubstructMatch(p_menthene_pattern1) or
            mol.HasSubstructMatch(p_menthene_pattern2)):
        return True, "Matches monoterpenoid criteria"

    # Additional check for rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    found_6_member = False
    found_5_member = False
    for ring in ring_info.AtomRings():
      if len(ring) == 6:
        found_6_member = True
      if len(ring) == 5:
         found_5_member = True

    if num_rings > 0 and not (found_6_member or found_5_member) :
         return False, "Monocyclic monoterpenoids usually contain a 5 or 6-membered ring"

    
    return False, "Does not match common monoterpenoid substructures or characteristics"