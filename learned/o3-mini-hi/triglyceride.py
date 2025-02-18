"""
Classifies: CHEBI:17855 triglyceride
"""
"""
Classifies: Triglycerides
Definition: Any glyceride resulting from the condensation of all three hydroxy groups of glycerol 
(propane-1,2,3-triol) with fatty acids.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride is defined as a glyceride in which all three hydroxy groups of glycerol
    are esterified with fatty acids.
    
    The following steps are used:
      1. Parse the SMILES string.
      2. Identify exactly 3 ester groups via the SMARTS pattern "[OX2][CX3](=[OX1])".
      3. For each ester group, retrieve the bridging oxygen atom and then find its neighbor
         (other than the carbonyl carbon) which is assumed to be part of the glycerol backbone.
      4. Check that exactly 3 unique glycerol carbons were identified.
      5. Verify that these three carbons form a connected chain:
         one central carbon (with 2 connections within the backbone) and two terminal carbons (1 connection each).
      6. Optionally, validate that the molecule has additional properties typical for triglycerides
         (e.g. molecular weight above 500 Da, sufficient carbon and oxygen counts and a high number of rotatable bonds).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a triglyceride, False otherwise.
        str: Reason for classification decision.
    """
    # Step 1: Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Step 2: Identify ester groups.
    # Ester groups are identified by the substructure -O-C(=O)- where
    # SMARTS "[OX2][CX3](=[OX1])" returns a tuple (O, C, O) for each match.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if ester_pattern is None:
        return False, "Error in ester SMARTS pattern."
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 3:
        return False, f"Found {len(ester_matches)} ester group(s); expected exactly 3 ester groups for a triglyceride."
    
    # Step 3: Identify the glycerol backbone atoms.
    # For each ester match, the atom at index 0 is the ester oxygen.
    # This oxygen will have two neighbors: one will be the carbonyl carbon (from the ester group match)
    # and the other is expected to be part of the glycerol backbone.
    glycerol_candidates = set()
    for match in ester_matches:
        # match[0] is the bridging oxygen; get the atom.
        oxygen_atom = mol.GetAtomWithIdx(match[0])
        # Get neighbors of the oxygen. One neighbor is already the ester carbon (match[1]).
        for neigh in oxygen_atom.GetNeighbors():
            # If this neighbor is a carbon but not the carbonyl carbon, consider it as glycerol candidate.
            if neigh.GetAtomicNum() == 6 and neigh.GetIdx() != match[1]:
                glycerol_candidates.add(neigh.GetIdx())
    
    if len(glycerol_candidates) != 3:
        return False, f"Expected 3 glycerol backbone carbons, found {len(glycerol_candidates)}."
    
    # Step 4: Check connectivity of the glycerol backbone.
    # In glycerol, the three carbons should be connected to one another in a chain:
    # two terminal carbons (each having 1 connection within the backbone)
    # and one central carbon (with 2 backbone connections).
    backbone_adj_counts = []
    # Convert set to list for iteration
    backbone_atoms = list(glycerol_candidates)
    for idx in backbone_atoms:
        atom = mol.GetAtomWithIdx(idx)
        # Count neighbors that are also in glycerol_candidates.
        count = 0
        for neigh in atom.GetNeighbors():
            if neigh.GetIdx() in glycerol_candidates:
                count += 1
        backbone_adj_counts.append(count)
    backbone_adj_counts.sort()
    if backbone_adj_counts != [1, 1, 2]:
        return False, "The three candidate glycerol carbons are not connected in a proper glycerol backbone (expected connectivity [1,1,2])."
    
    # Step 5: Additional property checks.
    # Check for a sufficiently high number of rotatable bonds (suggesting long fatty acid chains).
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Too few rotatable bonds, which suggests the fatty acid chains may be too short."
    
    # Check molecular weight (typically TG > 500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a triglyceride."
    
    # Check for a minimum count of oxygen and carbon atoms.
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, "Too few oxygen atoms for a triglyceride (expect at least 6)."
    
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbon atoms for a triglyceride (expect at least 20)."
    
    return True, "Molecule is classified as a triglyceride: it has a glycerol backbone with 3 esterified fatty acid chains."


# Example usage (uncomment to test):
# test_smiles = "O([C@H](COC(=O)CCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)C(=O)CCCCCCCCCCCC"
# result, message = is_triglyceride(test_smiles)
# print(result, message)