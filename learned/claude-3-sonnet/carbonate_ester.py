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
    # Carbonate ester pattern: R-O-C(=O)-O-R where R is carbon-based
    carbonate_pattern = Chem.MolFromSmarts("[C][O;X2][C;X3](=[O;X1])[O;X2][C]")
    
    # Monocarbonate pattern (R-O-C(=O)-OH)
    monocarbonate_pattern = Chem.MolFromSmarts("[C][O;X2][C;X3](=[O;X1])[O;X2][H]")
    
    # Cyclic carbonate pattern (specifically for 5-membered rings)
    cyclic_carbonate_pattern = Chem.MolFromSmarts("O=C1OC[C]1")
    
    # Exclusion patterns
    lactone_pattern = Chem.MolFromSmarts("O=C1[C,O][C,O][C,O]1")  # Cyclic esters
    anhydride_pattern = Chem.MolFromSmarts("[C;X3](=[O;X1])[O;X2][C;X3](=[O;X1])")
    carbamate_pattern = Chem.MolFromSmarts("[NX3][C;X3](=[O;X1])[O;X2]")
    ester_pattern = Chem.MolFromSmarts("[C;!$(C(=O)O)][C;X3](=[O;X1])[O;X2][C]")
    
    def is_valid_carbonate(match, mol):
        """Helper function to validate carbonate group"""
        c_atom = match[1]  # Central carbon
        c_atom_obj = mol.GetAtomWithIdx(c_atom)
        
        # Must have exactly 3 oxygens connected
        o_neighbors = sum(1 for atom in c_atom_obj.GetNeighbors() if atom.GetAtomicNum() == 8)
        if o_neighbors != 3:
            return False
            
        # Check that both O atoms are connected to carbon-based groups
        o1_atom = mol.GetAtomWithIdx(match[0])
        o2_atom = mol.GetAtomWithIdx(match[3])
        
        def has_carbon_neighbor(o_atom, exclude_idx):
            return any(n.GetAtomicNum() == 6 and n.GetIdx() != exclude_idx 
                      for n in o_atom.GetNeighbors())
        
        return (has_carbon_neighbor(o1_atom, c_atom) and 
                has_carbon_neighbor(o2_atom, c_atom))
    
    # Find matches
    carbonate_matches = mol.GetSubstructMatches(carbonate_pattern)
    monocarbonate_matches = mol.GetSubstructMatches(monocarbonate_pattern)
    
    # Get exclusion matches
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    anhydride_matches = mol.GetSubstructMatches(anhydride_pattern)
    carbamate_matches = mol.GetSubstructMatches(carbamate_pattern)
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Check linear carbonate matches
    for match in carbonate_matches:
        # Skip if part of excluded groups
        if any(match[1] in m for m in anhydride_matches):
            continue
        if any(match[1] in m for m in carbamate_matches):
            continue
        if any(match[1] in m for m in lactone_matches):
            continue
            
        if is_valid_carbonate(match, mol):
            return True, "Contains carbonate ester group (R-O-C(=O)-O-R)"
            
    # Check monocarbonate matches
    for match in monocarbonate_matches:
        if any(match[1] in m for m in anhydride_matches + lactone_matches):
            continue
        if is_valid_carbonate(match, mol):
            return True, "Contains monocarbonate ester group (R-O-C(=O)-OH)"
    
    return False, "No carbonate ester group found"