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
    carbon_check =  23 <= carbon_count <= 70
    
    # Check for isoprene units using a more general pattern
    isoprene_pattern = Chem.MolFromSmarts("[C]([C])([C])C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    isoprene_count = len(isoprene_matches)
    isoprene_check = isoprene_count >= 3

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    mol_wt_check = 300 <= mol_wt <= 1400
    
    # Check for branching (quaternary carbons with methyl groups)
    branching_pattern = Chem.MolFromSmarts("[C]([CH3])([CH3])([!H])")
    branching_matches = mol.GetSubstructMatches(branching_pattern)
    branching_check = len(branching_matches) >= 2

     # Check for a core substructure
    core_pattern = Chem.MolFromSmarts("C1CC(C)(C)CC2C(C)(C)CC3C(C)(C)CC(C)(C)C(C)(C)C(C)(C)C3CC2CC1")
    core_match = mol.HasSubstructMatch(core_pattern)

    # Combine criteria
    if carbon_check and mol_wt_check and branching_check and (isoprene_check or core_match):
        return True, "Matches criteria for a sesterterpenoid based on number of carbons, molecular weight, branching and isoprene units."
    else:
        reasons = []
        if not carbon_check:
             reasons.append(f"Number of carbons is {carbon_count}, expected between 23 and 70.")
        if not mol_wt_check:
            reasons.append(f"Molecular weight is {mol_wt}, expected between 300 and 1400.")
        if not branching_check:
            reasons.append(f"Too few branching points ({len(branching_matches)}), expected at least 2.")
        if not (isoprene_check or core_match):
             reasons.append(f"Too few isoprene units ({isoprene_count}) or core substructure was not found.")

        return False, "Does not match sesterterpenoid criteria: " + " ".join(reasons)