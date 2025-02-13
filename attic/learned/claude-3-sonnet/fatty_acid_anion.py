"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: CHEBI:29623 fatty acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the conjugate base of a fatty acid, with a deprotonated 
    carboxylic acid group (-COO-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylate group (-COO-)
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) != 1:
        return False, f"Must have exactly one carboxylate group, found {len(carboxylate_matches)}"

    # Count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 2:
        return False, "Carbon chain too short"
    
    # Check for multiple charges
    charge_sum = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if charge_sum != -1:
        return False, "Total charge must be -1"
        
    # Count other charged groups that would make this not a simple fatty acid anion
    other_anions = Chem.MolFromSmarts("[O-,N-,S-;!$([O-][C]=[O])]")
    if other_anions and mol.HasSubstructMatch(other_anions):
        return False, "Contains other anionic groups"
        
    # Look for multiple carboxylic/carboxylate groups
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,O-]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) > 1:
        return False, "Contains multiple carboxylic/carboxylate groups"

    # Check if the carboxylate is connected to carbon chain
    carboxylate_carbon = mol.GetAtomWithIdx(carboxylate_matches[0][0])
    neighbors = [n for n in carboxylate_carbon.GetNeighbors() 
                if n.GetAtomicNum() == 6]
    if not neighbors:
        return False, "Carboxylate not connected to carbon chain"

    # Optional: Check for common fatty acid features
    # Count rotatable bonds to verify chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2:
        return False, "Carbon chain too rigid/short for fatty acid"

    # Success case
    return True, "Contains single carboxylate group with appropriate carbon chain"