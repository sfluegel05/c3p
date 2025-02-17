"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
"""
Classifies: Aromatic primary alcohol
Definition:
  An aromatic primary alcohol is defined as any primary alcohol in which the alcoholic hydroxy group 
  is attached to a carbon which is itself bonded to an aromatic ring.
  
We improved the previous approach that solely relied on the SMARTS "c[CH2][OH]". Instead, we:
  1. Parse the molecule and add explicit hydrogens.
  2. Iterate over all carbon atoms.
  3. For each carbon, check if it is bonded to a terminal –OH group (i.e. its only oxygen neighbor has degree 1),
     and verify that it is a CH2 group (by counting its hydrogens).
  4. Also check whether that carbon is connected to at least one aromatic atom.
If all these criteria are met, we classify the molecule as an aromatic primary alcohol.
"""

from rdkit import Chem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    
    An aromatic primary alcohol is one in which the hydroxyl group (-OH) is attached to a primary
    carbon (CH2) and that carbon is also bonded to at least one aromatic atom. This implementation 
    explicitly checks for these criteria.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an aromatic primary alcohol, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens (this is crucial to properly count H atoms on the CH2 group)
    mol = Chem.AddHs(mol)
    
    # Iterate over all atoms to look for a primary alcohol group (CH2OH) whose 
    # carbon is bonded to at least one aromatic atom.
    for atom in mol.GetAtoms():
        # We are interested in carbon atoms only.
        if atom.GetAtomicNum() == 6:
            # Check if this carbon has an attached hydroxyl group and is CH2.
            oh_neighbor = None
            # Look at the neighbors.
            for nbr in atom.GetNeighbors():
                # Identify an oxygen that is part of a hydroxyl group.
                # Here we require that the oxygen is terminal (has degree 1, i.e. bonded only to the carbon)
                if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
                    oh_neighbor = nbr
                    break
            
            # If no hydroxyl neighbor was found, skip this carbon.
            if oh_neighbor is None:
                continue
            
            # Check that the carbon is primary: it should have exactly two hydrogen atoms.
            # Using GetTotalNumHs (which includes implicit hydrogens after AddHs)
            num_h = atom.GetTotalNumHs()
            if num_h != 2:
                continue
            
            # Check that the same carbon is bonded to at least one aromatic atom.
            has_aromatic_neighbor = False
            for nbr in atom.GetNeighbors():
                # Skip the oxygen we already used.
                if nbr.GetAtomicNum() == 8:
                    continue
                if nbr.GetIsAromatic():
                    has_aromatic_neighbor = True
                    break
            
            if has_aromatic_neighbor:
                return True, "Molecule contains an aromatic primary alcohol group"
    
    # If no matching substructure was found, return False.
    return False, "No aromatic primary alcohol group found"

# Example usage (for manual testing; can be removed if using as a module):
if __name__ == "__main__":
    examples = {
        "4-acetoxybenzyl alcohol": "CC(=O)OC1=CC=C(CO)C=C1",
        "3-pyridinemethanol": "C1=CC(=CN=C1)CO",
        "2-methylbenzyl alcohol": "CC1=C(CO)C=CC=C1",
        "4-cyanobenzyl alcohol": "C1(=CC=C(C=C1)C#N)CO",
        # Some additional examples to test borderline cases.
    }
    
    for name, smiles in examples.items():
        result, reason = is_aromatic_primary_alcohol(smiles)
        print(f"{name} ({smiles}): {result}, {reason}")