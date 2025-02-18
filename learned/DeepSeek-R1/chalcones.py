"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    Chalcones are 1,3-diphenylpropenones (Ar-CH=CH-CO-Ar) and derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Enhanced chalcone pattern: two aromatic rings connected via propenone (C=O-CH-CH)
    # Using more flexible SMARTS to account for different substituents and directions
    chalcone_pattern = Chem.MolFromSmarts("[c]C(=O)-C=C-[c]")
    matches = mol.GetSubstructMatches(chalcone_pattern)
    
    if not matches:
        # Try alternative pattern where double bond position is reversed (C=C-C=O)
        alt_pattern = Chem.MolFromSmarts("[c]-C=C-C(=O)-[c]")
        matches = mol.GetSubstructMatches(alt_pattern)
        if not matches:
            return False, "No chalcone core (Ar-C(=O)-C=C-Ar) found"

    # Verify the carbonyl group is a ketone (not ester/acid/amide)
    for match in matches:
        # Find carbonyl carbon index based on matched pattern
        # For [c]C(=O)-C=C-[c] pattern, carbonyl is at index 1
        # For [c]-C=C-C(=O)-[c] pattern, carbonyl is at index 3
        if len(match) == 5:  # First pattern has 5 atoms in SMARTS
            carbonyl_idx = match[1]
        else:  # Assume alt_pattern has 5 atoms
            carbonyl_idx = match[3]
        
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        carbon_neighbors = 0
        oxygen_found = False
        
        for neighbor in carbonyl_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors += 1
            elif neighbor.GetAtomicNum() == 8:
                # Verify oxygen is double-bonded and has no other connections
                bond = mol.GetBondBetweenAtoms(carbonyl_atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() != Chem.BondType.DOUBLE:
                    return False, "Carbonyl oxygen not properly bonded"
                if neighbor.GetDegree() > 1:
                    return False, "Carbonyl oxygen has additional bonds (ester/acid)"
                oxygen_found = True
            else:
                return False, "Non-carbon/oxygen atom attached to carbonyl"

        if not oxygen_found:
            return False, "Missing carbonyl oxygen"
        if carbon_neighbors != 2:
            return False, f"Carbonyl group has {carbon_neighbors} carbon neighbors (needs 2)"

    return True, "Contains chalcone core (Ar-C(=O)-C=C-Ar) with ketone group"