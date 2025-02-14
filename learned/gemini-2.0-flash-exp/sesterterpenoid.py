"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Descriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are C25 isoprenoids with possible rearrangements.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    carbon_check =  20 <= carbon_count <= 70
    
    # Check for isoprene units (C5H8)
    isoprene_pattern = Chem.MolFromSmarts("CC(C)=CC")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    isoprene_count = len(isoprene_matches)
    isoprene_check = isoprene_count >= 3

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    mol_wt_check = 300 <= mol_wt <= 900
    
    # Check number of rings
    ring_count = Descriptors.RingCount(mol)
    ring_check = ring_count >= 2

    #check for a core skeleton
    #core_pattern = Chem.MolFromSmarts("[C]1[C]([C])([C])[C]([C])([C])[C]2[C]1[C]([C])([C])[C]([C])([C])[C]3[C]2[C]([C])([C])[C]([C])([C])[C]3")
    #core_match = mol.HasSubstructMatch(core_pattern)

    # Combine criteria
    if carbon_check and isoprene_check and mol_wt_check and ring_check:
        return True, "Matches criteria for a sesterterpenoid based on number of carbons, isoprene units, molecular weight, and ring count."
    else:
        reasons = []
        if not carbon_check:
             reasons.append(f"Number of carbons is {carbon_count}, expected between 20 and 70.")
        if not isoprene_check:
             reasons.append(f"Too few isoprene units ({isoprene_count}). Expected at least 3.")
        if not mol_wt_check:
            reasons.append(f"Molecular weight is {mol_wt}, expected between 300 and 900.")
        if not ring_check:
            reasons.append(f"Ring count is {ring_count}, expected at least 2.")


        return False, "Does not match sesterterpenoid criteria: " + " ".join(reasons)