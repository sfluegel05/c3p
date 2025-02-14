"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdMolDescriptors

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a D-glucopyranose ring, notice the switched stereocenter at C2.
    glucose_pattern = Chem.MolFromSmarts("[C@@H]1([CH2]O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)O1")
    matches = mol.GetSubstructMatches(glucose_pattern)

    if not matches:
        return False, "No D-glucose ring found."
    
    for match in matches:
        # Get anomeric carbon (C1) index
        anomeric_carbon_index = match[0]
        anomeric_atom = mol.GetAtomWithIdx(anomeric_carbon_index)
        
        # Check neighbors for glycosidic bond (any non-ring atom)
        for neighbor in anomeric_atom.GetNeighbors():
            if neighbor.GetIdx() not in match: # Neighbor atom is NOT part of the ring
                #check if is hydrogen, do nothing if so
                if neighbor.GetAtomicNum() == 1:
                  continue 
                #check the chirality of the anomeric carbon. If it is S, then it is beta.
                if anomeric_atom.GetChiralTag() == Chem.ChiralType.CHI_S:
                    return True, "beta-D-glucoside detected (anomeric carbon has S configuration)."
    
    return False, "Not a beta-D-glucoside (no glycosidic bond with correct beta configuration)"