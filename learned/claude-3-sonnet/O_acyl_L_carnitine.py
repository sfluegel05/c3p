"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: O-acyl-L-carnitine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an O-acyl-L-carnitine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carnitine backbone with ester and carboxylate
    # [C@H,C@@H] captures both possible SMILES representations of L-carnitine
    carnitine_pattern = Chem.MolFromSmarts('[C@H,C@@H](CC([O-])=O)(CO[C,c,S])C[N+](C)(C)C')
    if not mol.HasSubstructMatch(carnitine_pattern):
        return False, "Missing L-carnitine backbone structure"
    
    # Check for ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts('[OX2][CX3](=[OX1])[#6]')
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "Missing ester linkage"
    
    # Check for trimethylammonium group
    trimethyl_pattern = Chem.MolFromSmarts('[NX4+](C)(C)(C)')
    if not mol.HasSubstructMatch(trimethyl_pattern):
        return False, "Missing trimethylammonium group"
    
    # Check for carboxylate group
    carboxylate_pattern = Chem.MolFromSmarts('C([O-])=O')
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "Missing carboxylate group"
    
    # Count number of charged groups (should have one + and one -)
    pos_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0)
    neg_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
    
    if pos_charge != 1 or neg_charge != -1:
        return False, f"Incorrect charge distribution: +{pos_charge}, {neg_charge}"
    
    # Additional check for deuterated versions
    deuterium_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[2H]')))
    if deuterium_count > 0:
        return True, "Valid O-acyl-L-carnitine (deuterated form)"
        
    return True, "Valid O-acyl-L-carnitine structure"