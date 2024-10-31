from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops

def is_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an amino acid zwitterion, defined as having a negatively 
    charged carboxyl group and a positively charged amino group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid zwitterion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for negatively charged carboxyl group [O-]C(=O)
    carboxyl_pattern = Chem.MolFromSmarts('[O-]C(=O)')
    has_carboxyl = len(mol.GetSubstructMatches(carboxyl_pattern)) > 0
    
    # Look for positively charged amino groups [NH3+] or [NH2+]
    amino_pattern1 = Chem.MolFromSmarts('[NH3+]')
    amino_pattern2 = Chem.MolFromSmarts('[NH2+]')
    has_amino = (len(mol.GetSubstructMatches(amino_pattern1)) > 0 or 
                len(mol.GetSubstructMatches(amino_pattern2)) > 0)

    if not has_carboxyl:
        return False, "Missing negatively charged carboxyl group"
    if not has_amino:
        return False, "Missing positively charged amino group"

    # Get formal charges
    pos_charges = 0
    neg_charges = 0
    for atom in mol.GetAtoms():
        charge = atom.GetFormalCharge() 
        if charge > 0:
            pos_charges += charge
        elif charge < 0:
            neg_charges += abs(charge)

    # Check if molecule has appropriate number of charges
    if pos_charges == 0 or neg_charges == 0:
        return False, "Missing required charges"
    
    return True, f"Contains negatively charged carboxyl group and positively charged amino group"
# Pr=1.0
# Recall=0.9811320754716981