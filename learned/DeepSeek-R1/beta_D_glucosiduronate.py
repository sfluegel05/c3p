"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: CHEBI: beta-D-glucosiduronate
"""
from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate has a beta-D-glucuronic acid moiety (deprotonated carboxylate)
    linked via an O-glycosidic bond, with correct stereochemistry and no additional substituents.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Define the beta-D-glucuronide substructure pattern with correct stereochemistry and glycosidic bond
    glucuronide_pattern = Chem.MolFromSmarts(
        "[O-]C(=O)[C@@H]1O[C@H]([C@H]([C@@H]([C@H]1O)O)O)O-*"
    )
    if glucuronide_pattern is None:
        return False, "Failed to parse glucuronide pattern"
    
    # Check for presence of the glucuronide group with glycosidic bond
    matches = mol.GetSubstructMatches(glucuronide_pattern)
    if not matches:
        return False, "No beta-D-glucuronide substructure found"
    
    # Check that the glucuronide group is connected via O-glycosidic bond (anomeric oxygen has external bond)
    for match in matches:
        # Get the anomeric oxygen (O in the glycosidic bond)
        # The pattern is [O-]C(=O)C@@H1O[C@H](...)
        # The oxygen after C@@H1 is the ring oxygen, and the next oxygen is the glycosidic oxygen
        # Assuming the glycosidic oxygen is the fourth atom in the match
        # Adjust indices based on SMARTS parsing
        # The SMARTS is [O-]C(=O)[C@@H]1O[C@H](...)
        # Atoms: 0: O-, 1: C(=O), 2: C@@H1, 3: O (ring), 4: C@H (anomeric carbon)
        # The glycosidic oxygen is atom 3 (ring oxygen) is part of the ring, but the anomeric carbon (4) is connected to O-R
        # Wait, need to verify the atom indices in the match
        # Alternatively, find the anomeric carbon and check its bonds
        anomeric_carbon = None
        for atom in glucuronide_pattern.GetAtoms():
            if atom.GetAtomMapNum() == 1:  # Using atom map numbers to track the anomeric carbon
                anomeric_carbon = atom.GetIdx()
                break
        if anomeric_carbon is None:
            continue
        
        # Get the corresponding atom in the molecule
        mol_anomeric_carbon = match[anomeric_carbon]
        # Check that the anomeric carbon has an oxygen connected to another atom (aglycone)
        neighbors = mol.GetAtomWithIdx(mol_anomeric_carbon).GetNeighbors()
        oxygen_neighbors = [n for n in neighbors if n.GetAtomicNum() == 8]
        if len(oxygen_neighbors) < 1:
            continue  # No oxygen connected to anomeric carbon
        # Check if any oxygen neighbor is connected to another atom outside the glucuronide
        for oxy in oxygen_neighbors:
            # Check if the oxygen is part of the glycosidic bond (connected to non-glucuronide atom)
            # Get all bonded atoms to this oxygen
            oxy_neighbors = oxy.GetNeighbors()
            # Exclude the anomeric carbon and check if there's another atom
            if any(n.GetIdx() not in match for n in oxy_neighbors):
                # Found a glycosidic bond
                return True, "Contains beta-D-glucuronide group with carboxylate and O-glycosidic bond"
    
    return False, "No valid beta-D-glucuronide substructure with O-glycosidic bond"