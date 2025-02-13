"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: CHEBI:61691 3-substituted propionyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-substituted propionyl-CoA(4-) based on its SMILES string.
    A 3-substituted propionyl-CoA(4-) is an acyl-CoA(4-) oxoanion with a 3-substituted propionyl group
    attached to the thiol group of coenzyme A, and a maximum acyl chain length of 20 carbons.

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

    # Look for CoA(4-) backbone with correct stereochemistry
    coa_pattern = Chem.MolFromSmarts("C(COP([O-])(=O)OP([O-])(=O)OC[C@@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12)(NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O))")
    coa_match = mol.GetSubstructMatch(coa_pattern)
    if not coa_match:
        return False, "Missing CoA(4-) backbone or incorrect stereochemistry"
        
    # Look for 3-substituted propionyl group attached via thioester
    propionyl_pattern = Chem.MolFromSmarts("C(=O)SC[C@@](C)(C)[CX4]")
    propionyl_matches = mol.GetSubstructMatches(propionyl_pattern)
    if len(propionyl_matches) != 1:
        return False, f"Found {len(propionyl_matches)} 3-substituted propionyl groups, need exactly 1"
    
    # Check for connectivity between CoA(4-) backbone and 3-substituted propionyl group
    coa_atoms = [mol.GetAtomWithIdx(idx) for idx in coa_match]
    propionyl_atom = mol.GetAtomWithIdx(propionyl_matches[0][2])
    coa_neighbors = [atom for atom in propionyl_atom.GetNeighbors() if atom in coa_atoms]
    if not coa_neighbors:
        return False, "CoA(4-) backbone and 3-substituted propionyl group are not connected"
    
    # Check for acyl chain length (maximum 20 carbons)
    acyl_chain = [propionyl_atom]
    current_atom = propionyl_atom
    while True:
        neighbors = [atom for atom in current_atom.GetNeighbors() if atom.GetAtomicNum() == 6]
        if len(neighbors) == 0:
            break
        elif len(neighbors) > 1:
            return False, "More than one path in the acyl chain"
        current_atom = neighbors[0]
        acyl_chain.append(current_atom)
    if len(acyl_chain) > 23:  # 20 carbons + carbonyl + alpha carbon + 3-substituted carbon
        return False, "Acyl chain is too long (> 20 carbons)"
    
    # Check for negative charge on phosphate groups
    negative_charges = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == -1)
    if negative_charges != 4:
        return False, "Incorrect number of negative charges (should be 4)"
    
    return True, "Contains CoA(4-) backbone with a 3-substituted propionyl group attached via thioester"