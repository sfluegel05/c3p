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

    # Basic carbon skeleton: Try to match isoprene-like units.
    isoprene_pattern = Chem.MolFromSmarts("[CX4]([CX4])([CX4])([CX4])~[CX4]")
    matches_isoprene = mol.GetSubstructMatches(isoprene_pattern)

    # Check molecular weight: A monoterpene (C10H16) has a MW around 136
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 350: # Adjust upper bound, but keep the maximum molecular weight under 350
       return False, "Molecular weight too high for a monoterpenoid."

    # Check carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8 or c_count > 12 : #Allow for loss of some methyl groups
        return False, f"Too few or too many carbons: {c_count}. Monoterpenoids typically have around 10 carbons."

    # Check for presence of oxygen containing functional groups (alcohol, ketone, ether or ester)
    oxygen_pattern = Chem.MolFromSmarts("[OX2,OX1]")
    oxygen_matches = mol.GetSubstructMatches(oxygen_pattern)
    if len(oxygen_matches) == 0:
        return False, "Monoterpenoids typically contain oxygen functional groups."


    #Check ring counts: monoterpenoids can be acyclic, monocyclic, or bicyclic. We use the ring size information too.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 3:
        return False, "Monoterpenoids rarely contain more than 2-3 rings"

    # Check for 6-membered rings
    found_6_member = False
    for ring in ring_info.AtomRings():
      if len(ring) == 6:
        found_6_member = True
        break

    # If there is no isoprene match but there are one or more rings, then a 6-membered ring is required.
    if (len(matches_isoprene) == 0 and num_rings > 0 and not found_6_member):
        return False, "Monocyclic and bicyclic monoterpenoids usually include at least one six membered ring"
    
    return True, "Matches monoterpenoid criteria"