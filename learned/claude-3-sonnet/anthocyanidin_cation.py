"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: CHEBI:51744 anthocyanidin cation
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is an oxygenated derivative of flavylium (2-phenylchromenylium).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for positive charge
    if sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == 1) != 1:
        return False, "Molecule does not have a single positive charge"

    # Check for presence of rings (anthocyanidins are polycyclic)
    if Chem.GetSSSR(mol) == []:
        return False, "Molecule does not contain rings"

    # Check for presence of oxygen atoms
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 3:
        return False, "Molecule does not contain enough oxygen atoms"

    # Check for presence of aromatic rings
    aromatic_rings = [r for r in Chem.GetSymmSSSR(mol) if mol.GetRingInfo().IsAromaticRing(r)]
    if len(aromatic_rings) < 2:
        return False, "Molecule does not contain enough aromatic rings"

    # Check for presence of hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Molecule does not contain enough hydroxyl groups"

    # Check for molecular weight range (typically 200-500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, "Molecular weight outside the typical range for anthocyanidin cations"

    return True, "Molecule exhibits structural features consistent with an anthocyanidin cation"