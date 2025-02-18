"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: O-acyl-L-carnitine (CHEBI:74544)
"""
from rdkit import Chem
from rdkit.Chem import rdCIPLabeler

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    An O-acyl-L-carnitine has a carnitine backbone with L-configuration (R configuration) 
    and an O-acyl group attached via ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for core structure (without stereochemistry)
    # [N+] connected to any four groups (quaternary ammonium)
    # Connected to a carbon that has:
    # - Ester group (OC(=O)...)
    # - Carboxylate group (CC(=O)[O-])
    core_smarts = Chem.MolFromSmarts('[N+X4]-[C]-[C](OC(=O))-CC(=O)[O-]')
    matches = mol.GetSubstructMatches(core_smarts)
    
    if not matches:
        return False, "Core structure not found"
    
    # Check each candidate chiral center for R configuration
    for match in matches:
        # The chiral carbon is the third atom in the match (index 2)
        chiral_idx = match[2]
        atom = mol.GetAtomWithIdx(chiral_idx)
        
        # Must be a chiral center
        if atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED:
            continue
        
        # Compute CIP label
        try:
            rdCIPLabeler.AssignCIPLabels(mol)
            cip_label = atom.GetProp('_CIPCode')
        except:
            continue
        
        if cip_label == 'R':
            return True, "Core structure with R configuration (L-carnitine) found"
    
    return False, "No core structure with R configuration found"