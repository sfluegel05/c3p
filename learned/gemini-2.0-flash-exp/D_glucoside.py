"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside is a molecule with a D-glucose moiety linked to another molecule via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
         bool: True if the molecule is a D-glucoside, False otherwise.
         str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Generalized SMARTS pattern for D-glucose (covers both linear and cyclic forms), does not include the glycosidic oxygen
    glucose_smarts = "[CX4]([OX2])([CX4])([CX4])[CX4]([OX2])[CX4]([OX2])"
    glucose_pattern = Chem.MolFromSmarts(glucose_smarts)
    
    # Find all substructure matches of D-glucose
    glucose_matches = mol.GetSubstructMatches(glucose_pattern)
    
    if not glucose_matches:
       return False, "No glucose moiety found"

    for match in glucose_matches:
        for anomeric_carbon_idx in match:
             anomeric_carbon = mol.GetAtomWithIdx(anomeric_carbon_idx)
             
             for neighbor in anomeric_carbon.GetNeighbors():
                 if neighbor.GetAtomicNum() == 8: # Check for Oxygen
                     
                    # Check if the Oxygen is linked to a non-glucose atom
                    is_glycosidic = False
                    for neighbor_of_neighbor in neighbor.GetNeighbors():
                         if neighbor_of_neighbor.GetIdx() not in match:
                            is_glycosidic = True
                            glycosidic_carbon_idx = neighbor_of_neighbor.GetIdx()
                            break
                    
                    if is_glycosidic:
                        # Check stereochemistry of the anomeric carbon
                        
                        # Get the anomeric C in the pattern
                        anomeric_carbon_in_pattern = 0
                        # Get the chiral carbon in the D-glucose pattern
                        chiral_carbon_in_pattern = 1
                        
                        anomeric_carbon_match = mol.GetAtomWithIdx(anomeric_carbon_idx)
                        chiral_carbon_match = mol.GetAtomWithIdx(match[chiral_carbon_in_pattern])

                        # The bond between the oxygen and the carbon of the other group
                        bond_to_other_group = mol.GetBondBetweenAtoms(neighbor.GetIdx(),glycosidic_carbon_idx)

                        if bond_to_other_group is None:
                           continue
                         
                        # Get neighbors of the anomeric carbon. This must be the non-glycosidic bond part
                        anomeric_neighbors = [ n.GetIdx() for n in anomeric_carbon_match.GetNeighbors() if n.GetIdx() != neighbor.GetIdx()]
                        if len(anomeric_neighbors) != 3:
                           continue
                        
                        # Find the highest priority atom connected to the anomeric carbon
                        highest_priority_atom_idx = -1
                        highest_priority = -1
                        for neighbor_idx in anomeric_neighbors:
                            priority = mol.GetAtomWithIdx(neighbor_idx).GetAtomicNum()
                            if priority > highest_priority:
                                highest_priority = priority
                                highest_priority_atom_idx = neighbor_idx

                        if highest_priority_atom_idx == -1:
                            continue # No heavy atom

                        #Find the hydrogen atom linked to the chiral carbon of the D-glucose
                        hydrogen_of_chiral_carbon = -1
                        for neighbor_idx in chiral_carbon_match.GetNeighbors():
                           if mol.GetAtomWithIdx(neighbor_idx).GetAtomicNum() == 1:
                              hydrogen_of_chiral_carbon = neighbor_idx
                              break
                        if hydrogen_of_chiral_carbon == -1:
                            continue

                        #Get the configuration of the chiral carbon of the glucopyranose in the glucose unit
                        configuration_of_chiral_carbon = anomeric_carbon_match.GetChiralTag()

                        # Check if the stereochemistry is correct for D-glucose.
                        # D-Glucose should have the configuration set to R at that chiral center
                        if configuration_of_chiral_carbon != Chem.ChiralType.CHI_TETRAHEDRAL_R:
                           continue
                        return True, "Contains a D-glucose moiety linked via a glycosidic bond."
    
    return False, "No D-glucose moiety linked via a glycosidic bond found"