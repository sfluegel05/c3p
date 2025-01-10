"""
Classifies: CHEBI:15705 L-alpha-amino acid
"""
"""
Classifies: L-alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

def is_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an L-alpha-amino acid based on its SMILES string.
    
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
        
    # Neutralize charges to standardize the structure
    uncharger = Chem.UnchargeModel()
    mol = uncharger.uncharge(mol)

    # Look for alpha-amino acid pattern: NH2-CH-COOH
    # Note: also matches charged forms NH3+ and COO-
    aa_pattern = Chem.MolFromSmarts('[NX3,NX4+][CX4H1][CX3](=[OX1])[OX1H0-,OX2H1]')
    matches = mol.GetSubstructMatches(aa_pattern)
    
    if not matches:
        return False, "No alpha-amino acid group found"
    
    # Check each potential alpha-amino acid center
    for match in matches:
        N_idx, C_idx, C_carb, O_idx = match
        
        # Get the alpha carbon atom
        alpha_carbon = mol.GetAtomWithIdx(C_idx)
        
        # Check if alpha carbon is chiral
        if alpha_carbon.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
            
        # Get neighbors of alpha carbon to check if it's truly chiral
        neighbors = [n.GetAtomicNum() for n in alpha_carbon.GetNeighbors()]
        if len(set(neighbors)) < 3:  # Need at least 3 different atoms for chirality
            continue
            
        # Check for L configuration (usually R in CIP rules due to priority)
        # We can use the presence of @@ in the SMILES to indicate R configuration
        atom_smiles = Chem.MolToSmiles(mol, rootedAtAtom=C_idx)
        
        # Most L-amino acids have R configuration due to CIP rules
        # (except cysteine and derivatives where -SH has higher priority)
        has_sulfur = 'S' in CalcMolFormula(mol)
        
        if ('@@' in atom_smiles and not has_sulfur) or ('@' in atom_smiles and has_sulfur):
            return True, "Found L-alpha-amino acid center with correct stereochemistry"
            
    return False, "No L-alpha-amino acid stereocenter found"