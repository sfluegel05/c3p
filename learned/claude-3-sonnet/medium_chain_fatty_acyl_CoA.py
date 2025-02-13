"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:36540 medium-chain fatty acyl-CoA
A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any medium-chain fatty acid.
"""

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Match Coenzyme A substructure components
    coenzyme_a_patterns = [
        Chem.MolFromSmarts("[N-]=[N+]=C1NC(=NC=N1)N"),  # Adenine base
        Chem.MolFromSmarts("OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12"),  # Ribose-phosphate backbone
        Chem.MolFromSmarts("C(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCS")  # Pantothenate arm
    ]
    if not all(mol.HasSubstructMatch(pat) for pat in coenzyme_a_patterns):
        return False, "Missing Coenzyme A substructure components"

    # Find acyl group attached to Coenzyme A
    acyl_group = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S' and atom.GetIsAromatic() and atom.GetDegree() == 2:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetAtomicNum() == 6 and neighbor.GetFormalCharge() == 0 and neighbor.GetIsAromatic():
                acyl_group = neighbor
                break

    if acyl_group is None:
        return False, "No acyl group found attached to Coenzyme A"

    # Traverse bonds to find fatty acid chain
    chain = []
    current = acyl_group
    while current.GetDegree() > 1:
        for neighbor in current.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic() == False and neighbor not in chain:
                chain.append(neighbor)
                current = neighbor
                break

    # Check chain length (6-12 carbons)
    chain_length = len(chain)
    if chain_length < 6 or chain_length > 12:
        return False, f"Fatty acid chain length is {chain_length}, expected 6-12"

    # Check for additional acyl groups, esters, alcohols
    additional_groups = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetFormalCharge() == 0 and atom.GetIsAromatic() == False and atom not in chain]
    if additional_groups:
        return False, "Additional acyl groups, esters, or alcohols present"

    # Check for unsaturation, substituents, heteroatoms in the chain
    unsaturated = any(atom.GetDegree() > 3 for atom in chain)
    substituted = any(atom.GetDegree() > 2 for atom in chain if atom.GetAtomicNum() == 6)
    heteroatoms = any(atom.GetAtomicNum() != 6 and atom.GetAtomicNum() != 1 for atom in chain)

    # Construct reason for classification
    reason = "Contains Coenzyme A substructure with a medium-chain fatty acyl group (6-12 carbons) attached."
    if unsaturated:
        reason += " The fatty acid chain is unsaturated."
    if substituted:
        reason += " The fatty acid chain is substituted."
    if heteroatoms:
        reason += " The fatty acid chain contains heteroatoms."

    return True, reason