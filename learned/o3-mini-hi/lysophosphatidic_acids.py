"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
#!/usr/bin/env python
"""
Classifies: lysophosphatidic acids

Definition:
    Any monoacylglycerol phosphate obtained by hydrolytic removal of one of the two acyl groups of any
    phosphatidic acid or derivatives therein.

Criteria used for classification:
    1. The molecule must contain a phosphate moiety.
    2. The molecule must contain a glycerol backbone.
    3. The molecule must exhibit exactly one ester bond typical for an acyl chain (monoacyl).
    
Note: This is a heuristic approach; due to structural variation some edge cases may not be fully captured.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines if a molecule is a lysophosphatidic acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as a lysophosphatidic acid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # 1. Check for the phosphate group.
    # This SMARTS looks for a phosphorus atom double-bonded to oxygen and bonded to two hydroxyl oxygen atoms.
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not found."
    
    # 2. Check for a glycerol backbone.
    # We use a simplified pattern for glycerol (three carbons with two terminal oxygen substituents and one internal OH):
    # This pattern ignores chiral tags and is meant as a rough pattern: O-C-C(O)-C-O.
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "Glycerol backbone not found."
    
    # 3. Identify the acyl ester moiety.
    # Lysophosphatidic acids are monoacyl, i.e. they have exactly one fatty acid attached via an ester bond.
    # The ester bond is represented by the substructure O-C(=O) (an oxygen connected to a carbonyl group).
    ester_pattern = Chem.MolFromSmarts("OC(=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups; expected exactly 1 for monoacyl glycerol phosphate."
    
    # Optionally, one might check for the presence of a sufficiently long fatty acid chain by counting carbon atoms 
    # in the substituent, or perform further validation on the overall molecular formula. As a heuristic,
    # we will check that the molecular weight reflects a lipid of expected size.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low to be a lysophosphatidic acid."
    
    return True, "Molecule contains a phosphate group, a glycerol backbone, and exactly one acyl ester group (monoacyl)."

# Example usage:
if __name__ == "__main__":
    # Test a sample lysophosphatidic acid SMILES from the provided examples.
    test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCC)(O)(O)=O"
    result, reason = is_lysophosphatidic_acids(test_smiles)
    print(result, reason)