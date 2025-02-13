"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: CHEBI:61691 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    A 3-substituted propionyl-CoA(4-) is an acyl-CoA(4-) oxoanion with a substituted propionyl group
    attached to the thiol group of coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-substituted propionyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for CoA(4-) backbone
    coa_pattern = Chem.MolFromSmarts("C(COP([O-])(=O)OP([O-])(=O)OC[C@@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12)(NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O))")
    coa_match = mol.GetSubstructMatch(coa_pattern)
    if not coa_match:
        return False, "Missing CoA(4-) backbone"
        
    # Look for substituted propionyl group attached via thioester
    propionyl_pattern = Chem.MolFromSmarts("C(=O)SC")
    propionyl_matches = mol.GetSubstructMatches(propionyl_pattern)
    if len(propionyl_matches) != 1:
        return False, f"Found {len(propionyl_matches)} propionyl groups, need exactly 1"
    propionyl_atom = mol.GetAtomWithIdx(propionyl_matches[0][0])
    
    # Check for substitution on propionyl group
    substitutions = sum(1 for atom in propionyl_atom.GetNeighbors() if atom.GetAtomicNum() != 8 and atom.GetAtomicNum() != 6)
    if substitutions < 1:
        return False, "Propionyl group is not substituted"
    
    # Check for connectivity between CoA(4-) backbone and substituted propionyl group
    coa_atoms = [mol.GetAtomWithIdx(idx) for idx in coa_match]
    coa_neighbors = [atom for atom in propionyl_atom.GetNeighbors() if atom in coa_atoms]
    if not coa_neighbors:
        return False, "CoA(4-) backbone and substituted propionyl group are not connected"
    
    # Check for negative charge on phosphate groups
    negative_charges = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == -1)
    if negative_charges != 4:
        return False, "Incorrect number of negative charges (should be 4)"
    
    return True, "Contains CoA(4-) backbone with a substituted propionyl group attached via thioester"