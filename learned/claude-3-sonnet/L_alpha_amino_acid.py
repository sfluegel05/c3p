"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: CHEBI:33709 L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    An L-alpha-amino acid has L-configuration at the alpha-carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an L-alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for alpha-amino acid pattern: N[C@H](*)C(=O)O
    amino_acid_pattern = Chem.MolFromSmarts("N[C@H](*)[C@@](=O)O")
    matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    if not matches:
        return False, "No alpha-amino acid substructure found"
    
    # Get the alpha-carbon atom
    alpha_carbon_idx = matches[0][1]
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
    
    # Check for L-configuration at alpha-carbon
    if alpha_carbon.GetProp("_ChiralityPossible") == "No":
        return False, "Alpha-carbon is not a chiral center"
    
    # Get atom indices of atoms attached to alpha-carbon
    atom_map = {}
    for bond in alpha_carbon.GetBonds():
        idx = bond.GetOtherAtomIdx(alpha_carbon_idx)
        atom_map[bond.GetStereo()] = idx
    
    # Determine stereochemistry based on bond directions
    if atom_map.get(Chem.BondStereo.BEGINWEDGE) < atom_map.get(Chem.BondStereo.BEGINDASH):
        return True, "Contains an L-alpha-amino acid substructure"
    else:
        return False, "Contains a D-alpha-amino acid substructure"