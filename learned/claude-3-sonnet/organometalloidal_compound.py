"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    An organometalloidal compound has bonds between metalloid atoms (primarily As, but also B, Si, Ge, Sb, Te)
    and carbon atoms of an organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find arsenic atoms
    as_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 33]
    if not as_atoms:
        return False, "No arsenic atoms found"

    # Helper function to check if carbon is part of an organyl group
    def is_organyl_carbon(carbon_atom):
        # Check if it's part of common organic groups
        organic_pattern = True
        
        # Not organyl if carbon is part of carbonate or similar
        oxygen_neighbors = sum(1 for n in carbon_atom.GetNeighbors() 
                             if n.GetAtomicNum() == 8)
        if oxygen_neighbors >= 3:
            organic_pattern = False
            
        # Check for reasonable valence
        if carbon_atom.GetTotalValence() > 4:
            organic_pattern = False
            
        return organic_pattern

    valid_as_c_bonds = 0
    
    for as_atom in as_atoms:
        # Get formal charge and oxidation state indicators
        formal_charge = as_atom.GetFormalCharge()
        
        # Get carbon neighbors
        carbon_neighbors = [n for n in as_atom.GetNeighbors() 
                          if n.GetAtomicNum() == 6]
        
        # Check each carbon neighbor
        for carbon in carbon_neighbors:
            if is_organyl_carbon(carbon):
                # Verify bond type
                bond = mol.GetBondBetweenAtoms(as_atom.GetIdx(), carbon.GetIdx())
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    valid_as_c_bonds += 1
                    
    # Additional validation checks
    if valid_as_c_bonds == 0:
        return False, "No valid As-C bonds to organyl groups found"
        
    # Check for problematic patterns that might indicate inorganic species
    problematic_patterns = [
        "[As+]([O-])([O-])([O-])",  # arsenate
        "[As]([O-])([O-])([O-])",   # arsenite
        "[As](=O)([O-])([O-])",     # arsenate
    ]
    
    for pattern in problematic_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            return False, "Contains inorganic arsenic pattern"
            
    # Check total number of bonds to arsenic
    for as_atom in as_atoms:
        if as_atom.GetTotalNumHs() + len(as_atom.GetBonds()) > 5:
            return False, "Invalid arsenic valence"
            
    # Success case
    return True, f"Contains {valid_as_c_bonds} valid As-C bond(s) to organyl groups"