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

    # Define patterns for different types of carbonate esters
    
    # Linear carbonate pattern (R-O-C(=O)-O-R)
    linear_carbonate = Chem.MolFromSmarts("[#6]-[#8X2]-[#6X3](=[#8X1])-[#8X2]-[#6]")
    
    # Monocarbonate pattern (R-O-C(=O)-OH)
    monocarbonate = Chem.MolFromSmarts("[#6]-[#8X2]-[#6X3](=[#8X1])-[#8X2]-[H]")
    
    # Cyclic carbonate patterns
    # General cyclic carbonate
    cyclic_carbonate = Chem.MolFromSmarts("O=C1O[C,N][C,O,N]O1")
    
    # 5-membered cyclic carbonate
    cyclic_5 = Chem.MolFromSmarts("O=C1OC[C,O]O1")
    
    # 6-membered cyclic carbonate
    cyclic_6 = Chem.MolFromSmarts("O=C1OC[C,O][C,O]O1")
    
    # Fused cyclic carbonate
    fused_carbonate = Chem.MolFromSmarts("O=C1O[C]2[C,O][C,O][C,O][C,O]2O1")
    
    # Exclusion patterns
    carbamate = Chem.MolFromSmarts("[NX3][C;X3](=[O;X1])[O;X2]")
    anhydride = Chem.MolFromSmarts("[CX3](=O)[OX2][CX3](=O)")
    
    def check_carbonate_center(mol, match):
        """Verify that the carbonate carbon has correct connectivity"""
        c_idx = match[2]  # Central carbon index (adjust based on pattern)
        c_atom = mol.GetAtomWithIdx(c_idx)
        
        # Must have exactly 3 oxygens connected
        o_neighbors = sum(1 for atom in c_atom.GetNeighbors() if atom.GetAtomicNum() == 8)
        return o_neighbors == 3

    # Check for matches
    if mol.HasSubstructMatch(linear_carbonate):
        matches = mol.GetSubstructMatches(linear_carbonate)
        for match in matches:
            if check_carbonate_center(mol, match):
                return True, "Contains linear carbonate ester group (R-O-C(=O)-O-R)"
                
    if mol.HasSubstructMatch(monocarbonate):
        matches = mol.GetSubstructMatches(monocarbonate)
        for match in matches:
            if check_carbonate_center(mol, match):
                return True, "Contains monocarbonate ester group (R-O-C(=O)-OH)"
    
    if mol.HasSubstructMatch(cyclic_carbonate):
        # Check it's not a carbamate or anhydride
        if not (mol.HasSubstructMatch(carbamate) or mol.HasSubstructMatch(anhydride)):
            return True, "Contains cyclic carbonate ester group"
            
    if mol.HasSubstructMatch(cyclic_5):
        return True, "Contains 5-membered cyclic carbonate ester"
        
    if mol.HasSubstructMatch(cyclic_6):
        return True, "Contains 6-membered cyclic carbonate ester"
        
    if mol.HasSubstructMatch(fused_carbonate):
        return True, "Contains fused cyclic carbonate ester"

    return False, "No carbonate ester group found"