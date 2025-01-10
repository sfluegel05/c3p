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
    An O-acyl-L-carnitine has an L-carnitine backbone with an acyl group attached via an ester bond.
    
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

    # Check for basic carnitine backbone structure:
    # - A carbon with an ester group
    # - Connected to a CH2 with carboxylate
    # - Connected to a CH2 with trimethylammonium
    # Note: [C@H] or [C@@H] represents the L-configuration
    carnitine_pattern = Chem.MolFromSmarts('[C@H,C@@H](CC([O-])=O)(CO)C[N+](C)(C)C')
    
    if not mol.HasSubstructMatch(carnitine_pattern):
        return False, "Missing L-carnitine backbone structure"
        
    # Verify the ester linkage is present
    # This pattern matches any acyl group (-C(=O)R) attached to oxygen
    ester_pattern = Chem.MolFromSmarts('OC(=O)[#6]')
    if not mol.HasSubstructMatches(ester_pattern):
        return False, "Missing ester linkage"

    # Check for trimethylammonium group
    # Allow for deuterated methyl groups
    trimethyl_pattern = Chem.MolFromSmarts('[N+]([CH3,CD3])([CH3,CD3])([CH3,CD3])')
    if not mol.HasSubstructMatch(trimethyl_pattern):
        return False, "Missing trimethylammonium group"
    
    # Check for carboxylate group
    carboxylate_pattern = Chem.MolFromSmarts('CC([O-])=O')
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "Missing carboxylate group"
    
    # Verify charges
    pos_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0)
    neg_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0)
    
    if pos_charge != 1 or neg_charge != -1:
        return False, f"Incorrect charge distribution: +{pos_charge}, {neg_charge}"

    # Check stereochemistry
    # In L-carnitine, the ester group and the trimethylammonium-methylene group 
    # should be on opposite sides of the central carbon
    chiral_centers = Chem.FindMolChiralCenters(mol)
    if not chiral_centers:
        return False, "Missing required stereocenter"
    
    # Additional check for deuterated versions
    deuterium_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[2H]')))
    if deuterium_count > 0:
        return True, "Valid O-acyl-L-carnitine (deuterated form)"
        
    return True, "Valid O-acyl-L-carnitine structure"