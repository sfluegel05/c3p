"""
Classifies: CHEBI:28874 phosphatidylinositol
"""
"""
Classifies: Phosphatidylinositol
Definition: Any glycerophosphoinositol having one phosphatidyl group 
esterified to one of the hydroxy groups of inositol.
We look for:
  1. An inositol ring (for example, the typical myo-inositol substructure).
  2. A P atom (phosphorus) that is connected via oxygen to that ring.
  3. The presence of at least two ester groups (OC(=O)) that would indicate
     the fatty acyl chains of the phosphatidyl group.
  4. A rough molecular weight cutoff to avoid very small molecules.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylinositol(smiles: str):
    """
    Determines if a molecule is a phosphatidylinositol based on its SMILES string.
    
    A phosphatidylinositol is characterized by a myo-inositol ring that is
    phosphorylated, with the phosphate ester further linked to a diacylglycerol
    (phosphatidyl) group. Thus, this function first checks the presence of an
    inositol ring, then looks for a phosphorus atom having an oxygen connection
    to the inositol and at least two ester groups (OC(=O)) that suggest fatty acyl chains.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as a phosphatidylinositol, otherwise False.
        str: A reason summarizing the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens (helpful for substructure matching of -OH groups)
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a myo-inositol ring.
    # This pattern covers a typical myo-inositol: a six-membered ring with 5 hydroxyl groups.
    inositol_smarts = "C1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)1"
    inositol_pattern = Chem.MolFromSmarts(inositol_smarts)
    inositol_matches = mol.GetSubstructMatches(inositol_pattern)
    if not inositol_matches:
        return False, "No inositol ring found"
    
    # Define a SMARTS for a fatty acyl ester group as OC(=O)
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, "Fewer than 2 ester (OC(=O)) groups found; insufficient for a phosphatidyl group"

    # Examine phosphorus atoms that may serve as the phosphoryl bridge.
    # We want a phosphorus that is linked to an oxygen coming from the inositol ring.
    phosphate_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # Phosphorus
            neighbors = atom.GetNeighbors()
            # Check if any neighbor is oxygen and is part of an inositol match
            inositol_attachment = False
            for nbr in neighbors:
                if nbr.GetAtomicNum() == 8:  # Oxygen
                    # Check if the oxygen atom is in any of the inositol substructure matches
                    for match in inositol_matches:
                        if nbr.GetIdx() in match:
                            inositol_attachment = True
                            break
                    if inositol_attachment:
                        break
            if not inositol_attachment:
                continue  # Try next P atom
            
            # For a proper phosphatidyl group, the phosphorus typically carries
            # several oxygens (usually 3 or more substituents) and is connected to the glycerol arm.
            oxy_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 8]
            if len(oxy_neighbors) < 3:
                continue  # Not enough oxygen substituents
            
            # Optional: Check molecular weight (phosphatidylinositols are relatively large lipids)
            mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
            if mol_wt < 600:
                continue  # Too low in molecular weight for a genuine PI
            
            # If all conditions are met, consider the molecule as phosphatidylinositol.
            phosphate_found = True
            return True, "Contains myo-inositol ring linked via phosphate to a diacylglycerol (phosphatidyl) group"
    
    if not phosphate_found:
        return False, "No appropriate phosphate bridge linking an inositol ring to a phosphatidyl group was found"
    
    # Fallback (should not normally be reached)
    return False, "Molecule does not meet the criteria for phosphatidylinositol"
    
# Example usage (uncomment for testing):
# test_smiles = "P(OC1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)1)(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCC)(O)=O"
# print(is_phosphatidylinositol(test_smiles))