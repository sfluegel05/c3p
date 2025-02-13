"""
Classifies: CHEBI:46722 carbonate ester
"""
"""
Classifies: carbonate ester
Definition: Any carbonate that is carbonic acid in which the hydrogens have been replaced by organyl groups
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_carbonate_ester(smiles: str):
    """
    Determines if a molecule is a carbonate ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a carbonate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define patterns
    # Linear carbonate ester pattern
    carbonate_pattern = Chem.MolFromSmarts("[O;X2][C;X3](=[O;X1])[O;X2]")
    
    # Cyclic carbonate pattern
    cyclic_carbonate_pattern = Chem.MolFromSmarts("O=C1OCC1") 
    
    # Monocarbonate ester pattern (R-O-C(=O)-OH)
    monocarbonate_pattern = Chem.MolFromSmarts("[O;X2][C;X3](=[O;X1])[O;H1]")
    
    # Exclusion patterns
    anhydride_pattern = Chem.MolFromSmarts("[C;X3](=[O;X1])[O;X2][C;X3](=[O;X1])")
    carbamate_pattern = Chem.MolFromSmarts("[NX3][C;X3](=[O;X1])[O;X2]")
    
    # Find matches
    carbonate_matches = mol.GetSubstructMatches(carbonate_pattern)
    cyclic_matches = mol.GetSubstructMatches(cyclic_carbonate_pattern)
    monocarbonate_matches = mol.GetSubstructMatches(monocarbonate_pattern)
    
    # Check exclusions
    anhydride_matches = mol.GetSubstructMatches(anhydride_pattern)
    carbamate_matches = mol.GetSubstructMatches(carbamate_pattern)
    
    if not (carbonate_matches or cyclic_matches or monocarbonate_matches):
        return False, "No carbonate ester group found"
    
    def is_valid_carbonate(match, mol):
        c_atom = match[1]  # Central carbon
        c_atom_obj = mol.GetAtomWithIdx(c_atom)
        
        # Must have exactly 3 oxygens connected
        o_neighbors = sum(1 for atom in c_atom_obj.GetNeighbors() if atom.GetAtomicNum() == 8)
        if o_neighbors != 3:
            return False
            
        # At least one oxygen must be connected to carbon (organic group)
        for o_idx in [match[0], match[2]]:
            o_atom = mol.GetAtomWithIdx(o_idx)
            for neighbor in o_atom.GetNeighbors():
                if (neighbor.GetAtomicNum() == 6 and 
                    neighbor.GetIdx() != c_atom and 
                    not any(neighbor.GetIdx() in m for m in anhydride_matches)):
                    return True
        return False
    
    # Check linear carbonate matches
    for match in carbonate_matches:
        # Skip if part of anhydride or carbamate
        if any(match[1] in m for m in anhydride_matches):
            continue
        if any(match[1] in m for m in carbamate_matches):
            continue
            
        if is_valid_carbonate(match, mol):
            return True, "Contains linear carbonate ester group (R-O-C(=O)-O-R)"
            
    # Check cyclic carbonate matches
    for match in cyclic_matches:
        c_atom = match[0]  # Central carbon
        if not any(c_atom in m for m in anhydride_matches):
            return True, "Contains cyclic carbonate ester group"
            
    # Check monocarbonate matches
    for match in monocarbonate_matches:
        if any(match[1] in m for m in anhydride_matches):
            continue
        if is_valid_carbonate(match, mol):
            return True, "Contains monocarbonate ester group (R-O-C(=O)-OH)"
    
    return False, "No valid carbonate ester group found"