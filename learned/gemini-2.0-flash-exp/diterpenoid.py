"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    Diterpenoids have a C20 skeleton (or modified) and are built from isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check carbon count (allow some loss of methyl groups from 20)
    if not 15 <= c_count <= 24:
        return False, f"Carbon count is {c_count}, not in the diterpenoid range (15-24)"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 450 :
         return False, f"Molecular weight of {mol_wt} is not in the diterpenoid range (250-450)"


    # Check for ring systems
    n_rings = rdMolDescriptors.CalcNumRings(mol)

    if n_rings < 1: #Most diterpenoids should have at least 2 rings.
         return False, f"Diterpenoids should have 2 or more rings, this one has {n_rings}"
    
    # Basic isoprene unit pattern - very general
    isoprene_pattern = Chem.MolFromSmarts("[CX4](C)([CX4])[CX4]~[CX4]") #Matches a branched isoprene unit

    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)

    if len(isoprene_matches) < 2:
        return False, "Too few isoprene units detected (less than 2)"
    
    
    
    return True, "Likely a diterpenoid based on carbon count, ring system, molecular weight and isoprene units"