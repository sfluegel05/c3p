"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of an amide group (-C(=O)N-)
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Expected exactly 1 amide group, found {len(amide_matches)}"
    
    # Get the carbonyl carbon of the amide group
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    
    if not carbonyl_matches:
        return False, "No carbonyl carbon found"

    carbonyl_atom_index = carbonyl_matches[0][0]
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_atom_index)
    
    # find all the carbon neighbours
    carbon_neighbors = [neighbor.GetIdx() for neighbor in carbonyl_atom.GetNeighbors() if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6]
    
    if not carbon_neighbors:
       return False, "No carbon attached to carbonyl"

    #Check for long chain attached to the carbonyl
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") #at least 6 carbons chain
    found_long_chain = False
    for neighbor_idx in carbon_neighbors:
        submol = Chem.PathToSubmol(mol, [carbonyl_atom_index, neighbor_idx])
        if submol.HasSubstructMatch(long_chain_pattern):
            #Check the length of the chain
            path = Chem.GetShortestPath(mol,carbonyl_atom_index, neighbor_idx)
            if len(path) > 6:
                found_long_chain = True
                break

    if not found_long_chain:
        return False, "No long carbon chain attached to the amide carbonyl"

    # Check for additional carbonyl directly on the nitrogen of the amide
    additional_carbonyl = Chem.MolFromSmarts("N[CX3](=[OX1])[!#1]") #check for additional carbonyl on N
    if mol.HasSubstructMatch(additional_carbonyl):
        return False, "Additional carbonyl group found on the amide nitrogen."

    # Check molecular weight and number of carbons for long chain requirement
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 4: # At least 4 carbon chain
        return False, "Not enough carbons for fatty chain"

    return True, "Molecule is a fatty amide"